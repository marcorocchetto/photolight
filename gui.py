
import sys

from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import numpy as np
import pickle
from scipy.stats import chi2

import matplotlib
from matplotlib.ticker import FuncFormatter,ScalarFormatter
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib as mpl
from matplotlib import rc
matplotlib.use('Qt4Agg',warn=False)
mpl.rcParams['axes.linewidth'] = 2 # set the value globally
mpl.rcParams['text.antialiased'] = True
rc('text', usetex=True) # use tex in plots
rc('font', **{'family':'serif','serif':['Palatino'],'size'   : 14})

sys.path.append('./classes')
sys.path.append('./library')

import parameters, target, dataset
from parameters import Parameters
from target import Target
from dataset import Dataset
from photometry import Photometry

import functions, img_scale, list_frames
from functions import *
from img_scale import *
from list_frames import *

gui_class = uic.loadUiType('gui.ui')[0]

plot_labels = ['None', 'Airmass', 'HJD', 'Signal to Noise',
               'X position', 'Y position', 'FWHM',
               'Ensemble flux', 'Ensemble flux total error',
               'Median subtracted magnitude',
               'Magnitude', 'Flux count',
               'Ensemble']

class Qt4MplCanvas(FigureCanvas):

    def __init__(self, parent):
        self.fig = Figure(facecolor='white')
        self.axes = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
        QtGui.QSizePolicy.Expanding,
        QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

class PlotWindow(QtGui.QMainWindow, gui_class):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle("Plot Window")
        self.main_widget = QtGui.QWidget(self)
        vbl = QtGui.QVBoxLayout(self.main_widget)
        self.qmc = Qt4MplCanvas(self.main_widget)
        ntb = NavigationToolbar(self.qmc, self.main_widget)
        vbl.addWidget(self.qmc)
        vbl.addWidget(ntb)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

class ApplicationWindow(QtGui.QMainWindow, gui_class):

    def __init__(self, parent=None):

        QtGui.QMainWindow.__init__(self)

        self.setupUi(self)

        # initialise plot window
        self.aw = PlotWindow()
        self.aw.show()


        self.p =self.load_photometry_object()

        # add apertures
        self.comboBox_apsize.addItems([str(x) for x in self.p.aperture.keys()])

        # add plot axis labels...
        self.comboBox_plotX.addItems(plot_labels)
        self.comboBox_plotY.addItems(plot_labels)
        self.comboBox_plotYerror.addItems(plot_labels)

        # frames tab
        self.load_treeWidget_frames()
        self.load_treeWidget_stars()



        self.plot()

        self.connect()

    def connect(self):
        self.pushButton_load_photometry.clicked.connect(self.load_photometry_object)

        self.doubleSpinBox_starn.valueChanged.connect(self.load_treeWidget_frames)
        self.treeWidget_frames.currentItemChanged.connect(self.plot)
        self.treeWidget_stars.currentItemChanged.connect(self.treeWidget_stars_currentItemChanged)
        self.treeWidget_frames.itemChanged.connect(self.treeWidget_frames_itemChanged)
        self.treeWidget_stars.itemChanged.connect(self.treeWidget_stars_itemChanged)
        self.comboBox_apsize.currentIndexChanged.connect(self.plot)
        self.comboBox_plotX.currentIndexChanged.connect(self.plot)
        self.comboBox_plotY.currentIndexChanged.connect(self.plot)
        self.comboBox_plotYerror.currentIndexChanged.connect(self.plot)
        self.doubleSpinBox_starn.valueChanged.connect(self.plot)
        self.checkBox_plot_excluded_frames.stateChanged.connect(self.plot)
        self.pushButton_badframes_run.clicked.connect(self.bad_frames_rejection)
        self.pushButton_optimize_alpha_obtain.clicked.connect(self.star_optimize_alpha_threshold)
        self.doubleSpinBox_optimize_run.clicked.connect(self.star_optimize_run)
        self.pushButton_saveplot_ascii.clicked.connect(self.save_plot_ascii)
        self.pushButton_saveplot_pdf.clicked.connect(self.save_plot_pdf)

        self.pushButton_errorbudget_plot.clicked.connect(self.plot_errorbudget)
        self.pushButton_errorbudget_saveascii.clicked.connect(self.save_errorbudget_ascii)
        self.pushButton_errorbudget_savepdf.clicked.connect(self.save_errorbudget_pdf)



    def star_optimize_alpha_threshold(self):
        alpha = self.doubleSpinBox_optimize_alpha.value()
        apsize = float(self.comboBox_apsize.currentText())
        nobs = np.count_nonzero(self.p.aperture[apsize].frames_mask)
        chi2dist = chi2(nobs-2)
        chi2limits = np.divide(chi2dist.interval(alpha), nobs-2)
        self.doubleSpinBox_optimize_lower.setValue(chi2limits[0])
        self.doubleSpinBox_optimize_upper.setValue(chi2limits[1])

    def star_optimize_run(self):
        apsize = float(self.comboBox_apsize.currentText())
        customlimits = [self.doubleSpinBox_optimize_lower.value(), self.doubleSpinBox_optimize_upper.value()]
        self.p.aperture[apsize].aplog = []
        rlm, star_exclude = self.p.aperture[apsize].optimize_ensemble(star_exclude=[], customlimits=customlimits)
        self.load_treeWidget_stars()
        self.plot()

    def bad_frames_rejection(self):
        starn = int(self.doubleSpinBox_starn.value())
        apsize = float(self.comboBox_apsize.currentText())
        outlier_threshold = [self.doubleSpinBox_badframes_lower.value(), self.doubleSpinBox_badframes_upper.value()]
        self.p.aperture[apsize].bad_frames_exclusion(starn=starn, outlier_threshold=outlier_threshold)
        print self.p.aperture[apsize].bad_frames_rlm_params

        self.load_treeWidget_frames()
        self.plot()

    def load_photometry_object(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, 'Select parameter file', 'Parfiles/')
        self.lineEdit_photometry_file.setText(filename)
        return pickle.load(open(filename))

    def load_treeWidget_frames(self):
        self.treeWidget_frames.clear()
        self.treeWidget_frames.setHeaderItem(QtGui.QTreeWidgetItem(['#','X', 'HJD', 'Normalised Flux', 'Flux error', 'SNR']))
        starn = int(self.doubleSpinBox_starn.value())
        apsize = float(self.comboBox_apsize.currentText())
        items = []
        for i in xrange(self.p.nobs):
            item_values = [str(i)]
            item_values.append(format_number(self.p.dataset.airmass[i], 3))
            item_values.append(format_number(self.p.dataset.hjd[i], 5))
            item_values.append(format_number(self.p.aperture[apsize].eflux[starn,i], 5))
            item_values.append(format_number(self.p.aperture[apsize].esigflux[starn,i], 5))
            item_values.append(format_number(self.p.aperture[apsize].signoise[starn,i], 2))
            item = QtGui.QTreeWidgetItem(item_values)
            if self.p.aperture[apsize].frames_mask[i]:
                item.setCheckState(0, QtCore.Qt.Checked)
            else:
                item.setCheckState(0, QtCore.Qt.Unchecked)
            items.append(item)
        self.treeWidget_frames.addTopLevelItems(items)

    def treeWidget_frames_itemChanged(self, item, column):
        idx = self.treeWidget_frames.indexOfTopLevelItem(item)
        apsize = float(self.comboBox_apsize.currentText())
        self.p.aperture[apsize].frames_mask[idx] = item.checkState(column)
        self.plot()

    def load_treeWidget_stars(self):
        self.treeWidget_stars.clear()
        self.treeWidget_stars.setHeaderItem(QtGui.QTreeWidgetItem(['#', 'Magnitude', 'Chi2']))
        apsize = float(self.comboBox_apsize.currentText())
        items = []
        for i in xrange(self.p.nstar):
            item_values = [str(i)]
            item_values.append(format_number(self.p.aperture[apsize].medmag[i], 5))
            item_values.append(format_number(self.p.aperture[apsize].rlm_redchi2[i], 4))
            item = QtGui.QTreeWidgetItem(item_values)
            if i in self.p.aperture[apsize].star_exclude:
                item.setCheckState(0, QtCore.Qt.Unchecked)
            else:
                item.setCheckState(0, QtCore.Qt.Checked)
            items.append(item)
        self.treeWidget_stars.addTopLevelItems(items)

    def treeWidget_stars_itemChanged(self, item, column):
        apsize = float(self.comboBox_apsize.currentText())
        root = self.treeWidget_stars.invisibleRootItem()
        star_exclude = []
        for i in range(root.childCount()):
            item = root.child(i)
            if item.checkState(0) == 0:
                star_exclude.append(i)
        self.p.aperture[apsize].magnitude_ensemble(star_exclude=star_exclude)
        self.load_treeWidget_stars()
        self.plot()

    def treeWidget_stars_currentItemChanged(self, item, column):
        self.doubleSpinBox_starn.setValue(self.treeWidget_stars.currentIndex().row())

    def save_plot_pdf(self):

        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save plot', '')
        self.aw.qmc.fig.savefig(str(filename), dpi=80,  bbox_inches='tight')

    def save_plot_ascii(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save to ascii', '')

        apsize = float(self.comboBox_apsize.currentText())
        x = self.get_plot_values(self.comboBox_plotX.currentIndex())
        y = self.get_plot_values(self.comboBox_plotY.currentIndex())
        yerr = self.get_plot_values(self.comboBox_plotYerror.currentIndex())
        mask = self.p.aperture[apsize].frames_mask
        if isinstance(x, (np.ndarray, np.generic)) and isinstance(y, (np.ndarray, np.generic)):

            if not isinstance(yerr, (np.ndarray, np.generic)):
                yerr = np.zeros(len(x))

            out = np.zeros((len(x),3))
            out[:,0] = x
            out[:,1] = y
            out[:,2] = yerr
            np.savetxt('%s_full.dat' % filename, out)
            out = np.zeros((len(x[~mask]),3))
            out[:,0] = x[~mask]
            out[:,1] = y[~mask]
            out[:,2] = yerr[~mask]
            np.savetxt('%s_excluded.dat' % filename, out)
            out = np.zeros((len(x[mask]),3))
            out[:,0] = x[mask]
            out[:,1] = y[mask]
            out[:,2] = yerr[mask]
            np.savetxt('%s_clean.dat' % filename, out)

    def plot(self, ploterror=False):
        apsize = float(self.comboBox_apsize.currentText())
        x = self.get_plot_values(self.comboBox_plotX.currentIndex())
        y = self.get_plot_values(self.comboBox_plotY.currentIndex())
        yerr = self.get_plot_values(self.comboBox_plotYerror.currentIndex())
        mask = self.p.aperture[apsize].frames_mask
        if isinstance(x, (np.ndarray, np.generic)) and isinstance(y, (np.ndarray, np.generic)):
            if not isinstance(yerr, (np.ndarray, np.generic)):
                yerr = np.zeros(len(x))
            self.aw.qmc.axes.clear()
            self.aw.qmc.axes.errorbar(x[mask], y[mask], yerr[mask], ls='none', fmt='o', color='blue')
            if self.checkBox_plot_excluded_frames.isChecked():
                self.aw.qmc.axes.errorbar(x[~mask], y[~mask], yerr[~mask], ls='none', fmt='o', color='red')
            idx = self.treeWidget_frames.currentIndex().row()
            if type(idx) is int:
                if idx > 0:
                    if isinstance(yerr, (np.ndarray, np.generic)):
                        self.aw.qmc.axes.errorbar(x[idx], y[idx], yerr[idx], ls='none', fmt='o', color='yellow', markersize=6)
                    else:
                        self.aw.qmc.axes.plot(x[idx], y[idx], 'o', color='yellow', markersize=6)
            self.aw.qmc.axes.set_xlabel(plot_labels[self.comboBox_plotX.currentIndex()])
            self.aw.qmc.axes.set_ylabel(plot_labels[self.comboBox_plotY.currentIndex()])
            set_backgroundcolor(self.aw.qmc.axes, 'white')
            self.aw.qmc.draw()
        if ploterror:
            try:
                self.aw2
                self.plot_errorbudget()
            except:
                pass

    def get_plot_values(self, idx):
        starn = int(self.doubleSpinBox_starn.value())
        apsize = float(self.comboBox_apsize.currentText())

        if plot_labels[idx] == 'Airmass':
            return self.p.dataset.airmass
        elif plot_labels[idx] == 'HJD':
            return self.p.dataset.hjd
        elif plot_labels[idx] == 'X position':
            return self.p.xc[apsize][:,starn]
        elif plot_labels[idx] == 'Y position':
            return self.p.yc[apsize][:,starn]
        elif plot_labels[idx] == 'FWHM':
            return self.p.dataset.hjd
        elif plot_labels[idx] == 'Signal to Noise':
            return self.p.aperture[apsize].signoise[starn]
        elif plot_labels[idx] == 'Ensemble flux':
            return self.p.aperture[apsize].eflux[starn]
        elif plot_labels[idx] == 'Ensemble flux total error':
            return self.p.aperture[apsize].esigflux[starn]
        elif plot_labels[idx] == 'Flux count':
            return self.p.fluxsum[apsize][:,starn]
        elif plot_labels[idx] == 'Magnitude':
            return self.p.aperture[apsize].mag[:,starn]
        elif plot_labels[idx] == 'Median subtracted magnitude':
            return self.p.aperture[apsize].med[:,starn]
        else:
            return 0

    def plot_errorbudget(self):
        try:
            self.aw2
        except:
            # initialise plot window
            self.aw2 = PlotWindow()
        apsize = float(self.comboBox_apsize.currentText())
        self.aw2.show()
        self.aw2.qmc.axes.clear()
        a = self.p.aperture[apsize]

        print a.stars_mask
        print a.star_exclude

        a.noise_sigscint_model = np.asarray(a.noise_sigscint_model)
        self.aw2.qmc.axes.scatter(a.medmag[a.stars_mask], a.noise_rms[a.stars_mask])
        self.aw2.qmc.axes.scatter(a.medmag[~a.stars_mask], a.noise_rms[~a.stars_mask], color='red')
        idx = self.treeWidget_stars.currentIndex().row()
        if type(idx) is int:
            if idx > 0:
                self.aw2.qmc.axes.plot(a.medmag[idx], a.noise_rms[idx], 'o', zorder=99, color='yellow')
        self.aw2.qmc.axes.plot(a.noise_mags_model, a.noise_sigmastar_model, label='Star noise')
        self.aw2.qmc.axes.plot(a.noise_mags_model, a.noise_sigmasky_model, label='Sky noise')
        self.aw2.qmc.axes.plot(a.noise_mags_model, a.noise_sigmaread_model, label='Readout noise')
        self.aw2.qmc.axes.plot(a.noise_mags_model, a.noise_sigscint_model, label='Scintillation noise')
        self.aw2.qmc.axes.plot(a.noise_mags_model, a.noise_sigscint_model/a.scintnoisecorr, label='Theoretical scintillation noise', ls='dashed')
        self.aw2.qmc.axes.plot(a.noise_mags_model,
                               np.sqrt(a.noise_sigmastar_model**2+a.noise_sigmasky_model**2+
                                       a.noise_sigmaread_model**2+a.noise_sigscint_model**2), label='Total noise')
        self.aw2.qmc.axes.legend(loc=2,frameon=False, prop={'size':12})
        self.aw2.qmc.axes.set_yscale("log", nonposy='clip')
        self.aw2.qmc.axes.set_xlabel('Instrumental magnitude')
        self.aw2.qmc.axes.set_ylabel('Root Mean Square (RMS)')
        set_backgroundcolor(self.aw2.qmc.axes, 'white')
        self.aw2.qmc.draw()

    def save_errorbudget_ascii(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, 'Save to ascii', '')
        apsize = float(self.comboBox_apsize.currentText())
        a = self.p.aperture[apsize]
        a.noise_sigscint_model = np.asarray(a.noise_sigscint_model)
        out = np.zeros((len(a.medmag[a.stars_mask]),2))
        out[:,0] = a.medmag[a.stars_mask]
        out[:,1] = a.noise_rms[a.stars_mask]
        np.savetxt('%s_rms_masked.dat' % filename, out)
        out = np.zeros((len(a.medmag[~a.stars_mask]),2))
        out[:,0] = a.medmag[~a.stars_mask]
        out[:,1] = a.noise_rms[~a.stars_mask]
        np.savetxt('%s_rms_unmasked.dat' % filename, out)
        out = np.zeros((len(a.medmag),2))
        out[:,0] = a.medmag
        out[:,1] = a.noise_rms
        np.savetxt('%s_rms.dat' % filename, out)
        out = np.zeros((len(a.noise_mags_model),2))
        out[:,0] = a.noise_mags_model
        out[:,1] = a.noise_sigmastar_model
        np.savetxt('%s_starnoise.dat' % filename, out)
        out[:,1] = a.noise_sigmasky_model
        np.savetxt('%s_skynoise.dat' % filename, out)
        out[:,1] = a.noise_sigmaread_model
        np.savetxt('%s_readoutnoise.dat' % filename, out)
        out[:,1] = a.noise_sigscint_model
        np.savetxt('%s_scintnoise_corr.dat' % filename, out)
        out[:,1] = a.noise_sigscint_model/a.scintnoisecorr
        np.savetxt('%s_scintnoise_not_corr.dat' % filename, out)
        out[:,1] = np.sqrt(a.noise_sigmastar_model**2+a.noise_sigmasky_model**2+
                                       a.noise_sigmaread_model**2+a.noise_sigscint_model**2)
        np.savetxt('%s_scintnoise_total_noise.dat' % filename, out)

    def save_errorbudget_pdf(self):
        try:
            self.aw2
            filename = QtGui.QFileDialog.getSaveFileName(self, 'Save error budget plot', '')
            self.aw2.qmc.fig.savefig(str(filename), dpi=80,  bbox_inches='tight')
        except:
            self.plot_errorbudget()
            self.save_errorbudget_pdf()

def format_number(value, dec):
    if isinstance(value, str):
        return value
    else:
        f = '{:.%if}' % dec
        return f.format(value)

def tick_formatter(x, p):

    if x < 1.0:
        return "%.1f" % x
    if x >= 1.0:
        return "%i" % x

# two useful functions to set background and foreground colours

def set_foregroundcolor(ax, color):
     '''For the specified axes, sets the color of the frame, major ticks,
         tick labels, axis labels, title and legend
     '''
     for tl in ax.get_xticklines() + ax.get_yticklines():
         tl.set_color(color)
     for spine in ax.spines:
         ax.spines[spine].set_edgecolor(color)
     for tick in ax.xaxis.get_major_ticks():
         tick.label1.set_color(color)
     for tick in ax.yaxis.get_major_ticks():
         tick.label1.set_color(color)
     ax.axes.xaxis.label.set_color(color)
     ax.axes.yaxis.label.set_color(color)
     ax.axes.xaxis.get_offset_text().set_color(color)
     ax.axes.yaxis.get_offset_text().set_color(color)
     ax.axes.title.set_color(color)
     lh = ax.get_legend()
     if lh != None:
         lh.get_title().set_color(color)
         lh.legendPatch.set_edgecolor('none')
         labels = lh.get_texts()
         for lab in labels:
             lab.set_color(color)
     for tl in ax.get_xticklabels():
         tl.set_color(color)
     for tl in ax.get_yticklabels():
         tl.set_color(color)


def set_backgroundcolor(ax, color):
     '''Sets the background color of the current axes (and legend).
         Use 'None' (with quotes) for transparent. To get transparent
         background on saved figures, use:
         pp.savefig("fig1.svg", transparent=True)
     '''
     ax.patch.set_facecolor(color)
     lh = ax.get_legend()
     if lh != None:
         lh.legendPatch.set_facecolor(color)





#put the main window here
def main():
    import sys
    qApp = QtGui.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()
    sys.exit(qApp.exec_())

if __name__ == '__main__':
    main()
