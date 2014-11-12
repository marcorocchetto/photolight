
class RobustStatistic():

    def __init__(self, x, y, sig=[], deg=1, chi2limit=0.9, customlimits=[]):

        from scipy.stats import chi2
        import numpy as np
        import statsmodels.api as sm
        import math

        # remove nan values. These might appear in any one of x, y or ysigma
        self.nanmask = np.logical_or(np.isnan(x), np.isnan(y))
        if len(sig) > 0: self.nanmask = np.logical_or(self.nanmask, np.isnan(sig))

        # get rid of nans
        x = x[~self.nanmask]
        y = y[~self.nanmask]
        if len(sig) > 0: sig = sig[~self.nanmask]

        # add intercept
        X = sm.add_constant(x)
        self.rlm = sm.RLM(y, X, missing='none', M=sm.robust.norms.TukeyBiweight()).fit()
        self.rlm_resid = self.rlm.resid  # residuals

        # compute mad statistics and get ratio
        self.mad = np.median(np.absolute(self.rlm_resid))
        self.sigmad = self.mad*1.4286
        self.ratio = self.rlm_resid/self.sigmad

        # mask values with ratio > 3.0 (reject bad outliers)
        self.mask = np.less(np.abs(self.ratio), 3.0)

        x = x[self.mask]
        y = y[self.mask]
        if len(sig):
            sig = sig[self.mask]

        # weight
        if len(sig) > 0:
            weight = 1/sig**2
        else:
            weight = None

        self.p = np.polyfit(x, y, deg, w=weight)

        # weighted linear fit to data - rejected points
        if len(sig) > 0:
            self.polyfit = np.polyfit(x, y, deg=1, w=1/sig**2)
        else:
            self.polyfit = np.polyfit(x, y, deg=1)

        self.afit = self.polyfit[0]
        self.bfit = self.polyfit[1]

        self.yfit = self.polyfit[1]+self.polyfit[0]*x
        self.x = x
        self.y = y
        if len(sig) > 0: self.sig = sig

        self.polyfit_resid = self.yfit - y
        self.polyfit_dof = len(x) - 2
        if len(sig) > 0:
            self.polyfit_zval = self.polyfit_resid/sig
            self.polyfit_chi2 = np.sum(self.polyfit_zval**2)
            self.polyfit_redchi2 = self.polyfit_chi2/self.polyfit_dof
        self.polyfit_rms = math.sqrt(np.sum(self.polyfit_resid**2)/self.polyfit_dof)
        self.polyfit_chi2dist = chi2(self.polyfit_dof)

        if len(customlimits) == 2:
            self.polyfit_chi2limits = customlimits
        else:
            self.polyfit_chi2limits = np.divide(self.polyfit_chi2dist.interval(chi2limit), float(self.polyfit_dof))

    def chi2_within_limits(self):
        if self.polyfit_redchi2 > self.polyfit_chi2limits[0] and self.polyfit_redchi2 < self.polyfit_chi2limits[1]:
            return True
        else:
            return False
