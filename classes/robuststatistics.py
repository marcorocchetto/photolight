from scipy.stats import chi2
import numpy as np
import statsmodels.api as sm
import math

class RobustStatistic(object):

    def __init__(self, x, y, sig=[], chi2limit=0.95, customlimits=[], outlier_threshold=[-3.0, 3.0], maxiter=50):

        nanmask = np.logical_or(np.isnan(x), np.isnan(y))

        mask = np.ones(len(x), dtype=bool)

        iter = 1
        rejn = 0

        # reject outliers

        while True:

            xrlm = x[mask & ~nanmask]
            yrlm = y[mask & ~nanmask]
            X = sm.add_constant(xrlm)
            rlm = sm.RLM(yrlm, X, missing='none', M=sm.robust.norms.TukeyBiweight()).fit()
            residuals = y - (rlm.params[0]+rlm.params[1]*x)
            mad = np.median(np.absolute(residuals))
            sigmad = mad*1.4286
            ratio = residuals/sigmad
            maskit = (ratio > outlier_threshold[0]) & (ratio < outlier_threshold[1])
            if np.array_equal(mask, maskit) == True or iter >= maxiter:
                self.outliers_mask = mask
                self.rlm_params = rlm.params
                self.niter = iter
                break
            else:
                mask = np.copy(maskit)

            iter += 1

        # weghted linear fit to cleaned data

        xlfit = x[self.outliers_mask]
        ylfit = y[self.outliers_mask]
        siglfit = sig[self.outliers_mask] if len(sig) else []
        weights = 1/siglfit**2 if len(siglfit) > 0 else None
        polyfit = np.polyfit(xlfit, ylfit, deg=1, w=weights)

        polyfit_resid = ylfit - (polyfit[1]+polyfit[0]*xlfit)
        polyfit_dof = len(xlfit) - 2

        if len(sig) > 0:
            polyfit_zval = polyfit_resid/siglfit
            polyfit_chi2 = np.sum(polyfit_zval**2)
            self.polyfit_redchi2 = polyfit_chi2/polyfit_dof

        self.polyfit_rms = math.sqrt(np.sum(polyfit_resid**2)/polyfit_dof)
        self.polyfit_chi2dist = chi2(polyfit_dof)
        self.polyfit = polyfit

        if len(customlimits) == 2:
            self.polyfit_chi2limits = customlimits
        else:
            self.polyfit_chi2limits = np.divide(self.polyfit_chi2dist.interval(chi2limit), float(polyfit_dof))

#
# class RobustStatistic():
#
#     def __init__(self, x, y, sig=[], deg=1, chi2limit=0.9, customlimits=[], outlier_threshold=3.0):
#
#         # remove nan values. These might appear in any one of x, y or ysigma
#         self.nanmask = np.logical_or(np.isnan(x), np.isnan(y))
#         if len(sig) > 0: self.nanmask = np.logical_or(self.nanmask, np.isnan(sig))
#
#         # get rid of nans
#         x = x[~self.nanmask]
#         y = y[~self.nanmask]
#         if len(sig) > 0: sig = sig[~self.nanmask]
#
#         # add intercept
#         X = sm.add_constant(x)
#         self.rlm = sm.RLM(y, X, missing='none', M=sm.robust.norms.TukeyBiweight()).fit()
#         self.rlm_resid = self.rlm.resid  # residuals
#
#         # compute mad statistics and get ratio
#         self.mad = np.median(np.absolute(self.rlm_resid))
#         self.sigmad = self.mad*1.4286
#         self.ratio = self.rlm_resid/self.sigmad
#
#         # mask values with ratio > 3.0 (reject bad outliers)
#         self.mask = np.less(np.abs(self.ratio), outlier_threshold)
#
#         x = x[self.mask]
#         y = y[self.mask]
#         if len(sig):
#             sig = sig[self.mask]
#
#         # weight
#         if len(sig) > 0:
#             weight = 1/sig**2
#         else:
#             weight = None
#
#         self.p = np.polyfit(x, y, deg, w=weight)
#
#         # weighted linear fit to data - rejected points
#         if len(sig) > 0:
#             self.polyfit = np.polyfit(x, y, deg=1, w=1/sig**2)
#         else:
#             self.polyfit = np.polyfit(x, y, deg=1)
#
#         self.afit = self.polyfit[0]
#         self.bfit = self.polyfit[1]
#         self.yfit = self.polyfit[1]+self.polyfit[0]*x
#         self.x = x
#         self.y = y
#         if len(sig) > 0: self.sig = sig
#
#         self.polyfit_resid = self.yfit - y
#         self.polyfit_dof = len(x) - 2
#         if len(sig) > 0:
#             self.polyfit_zval = self.polyfit_resid/sig
#             self.polyfit_chi2 = np.sum(self.polyfit_zval**2)
#             self.polyfit_redchi2 = self.polyfit_chi2/self.polyfit_dof
#         self.polyfit_rms = math.sqrt(np.sum(self.polyfit_resid**2)/self.polyfit_dof)
#         self.polyfit_chi2dist = chi2(self.polyfit_dof)
#
#         if len(customlimits) == 2:
#             self.polyfit_chi2limits = customlimits
#         else:
#             self.polyfit_chi2limits = np.divide(self.polyfit_chi2dist.interval(chi2limit), float(self.polyfit_dof))

    def chi2_within_limits(self):
        if self.polyfit_redchi2 > self.polyfit_chi2limits[0] and self.polyfit_redchi2 < self.polyfit_chi2limits[1]:
            return True
        else:
            return False
