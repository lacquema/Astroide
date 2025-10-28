#! /Users/lacquema/Astroide.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os
import numpy as np
from math import pi, sqrt
from random import random
import re
from scipy.stats import binned_statistic, binned_statistic_2d
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QCheckBox

# My packages
from WindowPlot import WindowPlot
from Parameters import *
from TransferData import TransferDataClass
from Resonance import ResClass
from Curve import CurveClass
from UtilsAnaSimu import find_delimiter

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d  # Import pour l'interpolation


### --- Tools Generating --- ###

class GeneralToolClass(QWidget):

    def __init__(self, ToolName, ToolStatus):
        super().__init__()

        # General parameters
        self.IndexSnap = 0 # Index of the snapshot
        self.colorList = ['black', 'blue', 'red', 'green', 'orange', 'pink']

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Layout
        self.Layout = QHBoxLayout()

        # Plot button
        self.BtnPlot = QPushButton(ToolName)
        self.BtnPlot.clicked.connect(self.Toggle_WindowPlot)
        self.BtnPlot.setStatusTip(ToolStatus)
        self.Layout.addWidget(self.BtnPlot)

        # Initialisation of plot windows
        self.WindowPlot = WindowPlot(ToolName)
        self.WindowPlot.SignalCloseWindowPlot.connect(lambda: self.BtnPlot.setEnabled(True))  # Enable button when window is closed

        # Connections between parameters and plot
        self.WindowPlot.WidgetParam.BtnReset.clicked.connect(self.ResetParams)
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.refresh_plots)

        # Widget container
        self.setLayout(self.Layout) # GeneralToolClass is directly the widget container

    # Open the plot window when the Plot button is clicked
    def Toggle_WindowPlot(self):
        """Open the plot window when the Plot button is clicked."""
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)
        self.refresh_plots()

    def refresh_plots(self):
        """Refresh all active plots when the refresh button is clicked."""
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            if WidgetPlot.isVisible():
                WidgetPlot.refresh_plot()

    def reset_plots(self):
        """Reset all active plots when the reset button is clicked."""
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            WidgetPlot.reset_plot()
        self.refresh_plots()
    
    # Change the index of the snapshot
    def Change_IndexSnap(self, value):
        self.IndexSnap = value
        self.refresh_plots()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        app.closeAllWindows()

    # Reset all widgets of the parameters window
    def ResetParams(self):
        """Reset all widgets of the parameters window."""
        for i in reversed(range(2, self.WindowPlot.WidgetParam.Layout.count())): 
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        self.InitParams()
        self.reset_plots()

    def Plot(self):
        return
    
    def InitParams(self):
        """Initialize parameters. This method should be overridden by subclasses."""
        return
    
    def replace_params_in_formula(self, formula, prefixe, nOrbitDefault):
        """Replace parameters and functions in the formula with their corresponding values."""
        # print(formula)
        # for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
        #     formula = formula.replace(f'[{num}]', f'[{str(int(num)-1)}]') # Replace [n] by [n-1]
        for param in ['t', 'a', 'e']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'{param}[{nOrbitDefault}]', formula) # Add [nOrbitDefault] to the parameter
            formula = re.sub(r'\b' + param + r'\b', f'{prefixe}{param}', formula) # Replace the parameter by its value
        for param in ['i', 'w', 'W']:
            formula = re.sub(r'\b' + param + r'\b(?!\[)', f'{param}[{nOrbitDefault}]', formula) # Add [nOrbitDefault] to the parameter
            formula = re.sub(rf'\b{param}\[(\d+)\]', rf'np.radians({prefixe}{param}[\1])', formula)
        for fonction in ['sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2', 'hypot', 'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh', 'exp', 'expm1', 'exp2', 'log', 'log10', 'log2', 'log1p', 'sqrt', 'square', 'cbrt', 'power', 'erf', 'erfc', 'gamma', 'lgamma', 'digamma', 'beta']:
            formula = re.sub(r'\b' + fonction + r'\b', f'np.{fonction}', formula) # Replace the function by its numpy equivalent
        # Convert the result to degrees if it is an angle in radians
        if any(f'np.{angle_func}' in formula for angle_func in ['arcsin', 'arccos', 'arctan', 'arctan2']):
            formula = f'np.degrees({formula})'
        # print('Formula is: '+formula)
        return formula
    
    def evaluate_formula(self, formula, prefix, nOrbitDefault=1):
        """Evaluate a formula."""
        print(formula)
        if not formula:
            return None
        formula = self.replace_params_in_formula(formula, prefix, nOrbitDefault)
        print(formula)
        try:
            return eval(formula)
        except Exception as e:
            print(f'Error evaluating formula: {e}')
            return None
        
    def UpdateParams(self):
        """Update parameters based on the current widget values."""
        # This method should be overridden by subclasses to update specific parameters
        pass

    def try_UpdateParams(self, WidgetPlot):
        """Try to update parameters and handle exceptions."""
        try:
            self.UpdateParams()
        except Exception as e:
            print('Wrong Parameters: ', e)
            WidgetPlot.Canvas.fig.draw()

    def LabelOf(self, var=str):
        """Return the label for a given variable."""
        labels = {
            't': 'Time',
            'P': 'Period',
            'a': 'Semi-major axis',
            'e': 'Eccentricity',
            'i': 'Inclinaison',
            'w': 'Argument of periastron',
            'W': 'Longitude of ascending node',
            'Chi2': 'Chi square',
            'M': 'Mean longitude'
        }
        return labels.get(var, 'Unknown variable')
    
    def UnitOf(self, var=str):
        """Return the unit for a given variable."""
        units = {
            't': '[Myr]',
            'P': '[yr]',
            'a': '[AU]',
            'e': '',
            'i': '[°]',
            'w': '[°]',
            'W': '[°]',
            'Chi2': '',
            'M': '[°]',
        }
        return units.get(var, 'Unknown variable')
    
    def convolve_positions(self, x, y, fwhm_au, extent):
        """
        Convolution 2D d'un nuage de particules en utilisant une PSF gaussienne spécifiée par sa FWHM.
        Le grid_size, extent et fwhm_per_pixel sont calculés automatiquement.

        Paramètres :
            x, y : array-like
                Coordonnées des particules (en AU)
            fwhm_au : float
                Largeur à mi-hauteur (FWHM) de la PSF (en AU)

        Retour :
            x_conv, y_conv : ndarray
                Coordonnées convoluées selon la PSF
        """
        # r = np.sqrt(x**2 + y**2)
        # r_max = np.max(r)
        # extent = np.ceil(r_max * (1 + 0.1))

        # Détermination automatique de fwhm_per_pixel
        if fwhm_au < 0.5:
            fwhm_per_pixel = 5
        elif fwhm_au < 3:
            fwhm_per_pixel = 4
        else:
            fwhm_per_pixel = 3

        grid_size = fwhm_au / fwhm_per_pixel
        sigma_pix = fwhm_au / (2.355 * grid_size)

        nbins = int(2 * extent / grid_size)
        x_edges = np.linspace(-extent, extent, nbins + 1)
        y_edges = np.linspace(-extent, extent, nbins + 1)

        density, _, _, _ = binned_statistic_2d(x, y, None, statistic='count', bins=[x_edges, y_edges])
        density = density.T

        smoothed_density = gaussian_filter(density, sigma=sigma_pix)

        x_centers = (x_edges[:-1] + x_edges[1:]) / 2
        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        flat_density = smoothed_density.ravel()
        flat_density /= flat_density.sum()

        N = len(x)
        idx = np.random.choice(len(flat_density), size=N, p=flat_density)
        x_idx = idx % smoothed_density.shape[1]
        y_idx = idx // smoothed_density.shape[1]
        x_conv = x_centers[x_idx] + np.random.uniform(-grid_size/2, grid_size/2, N)
        y_conv = y_centers[y_idx] + np.random.uniform(-grid_size/2, grid_size/2, N)

        return x_conv, y_conv


class SpaceView(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z):
        super().__init__('Space view', 'Space view of the system')

        # Data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.e_m = e_m
        self.Ex = Ex
        self.Ey = Ey
        self.Ez = Ez
        self.Epx = Epx
        self.Epy = Epy
        self.Epz = Epz
        self.X = X
        self.Y = Y
        self.Z = Z

        # Window plots initialisation
        self.WidgetPlotXY = self.WindowPlot.add_WidgetPlot(self.PlotXY, xlim=True, ylim=True)
        self.WidgetPlotXZ = self.WindowPlot.add_WidgetPlot(self.PlotXZ, xlim=True, ylim=True)
        self.WidgetPlotXYZ = self.WindowPlot.add_WidgetPlot(self.PlotXYZ, xlim=True, ylim=True, zlim=True, azim=True, elev=True)

        # Parameters initialisation
        self.InitParams()

    def InitParams(self):
        # Equal aspect ratio
        # self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter(Title='View :'))

        # Type of view
        self.ViewWidget = ComboBox('View', 'Dimension', ['face-on', 'edge-on', '3D'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ViewWidget)
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.ViewWidget.ComboParam.currentIndexChanged.connect(self.indexViewChanged)

        # Bodies' positions
        self.SizeBodies = 15
        self.CheckBodies = CheckBox("Bodies position")
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBodies)

        # self.SizeBodies = 15
        # self.SizeBodiesWidget = SpinBox('Size of bodies', 'Plot size of bodies', self.SizeBodies)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.SizeBodiesWidget)
        
        # self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked())
        # self.CheckBodies.CheckParam.stateChanged.connect(lambda: self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked()))

        # self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Orbits' tracing
        self.CheckOrbits = CheckBox("Bodies orbit")
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckOrbits)

        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter(Title='Representation :'))

        # Type of representation
        self.RepresWidget = ComboBox('Representation', 'Type of representation', ['scatter', 'density'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.RepresWidget)
        self.indexRepres = self.RepresWidget.ComboParam.currentIndex()
        self.RepresWidget.ComboParam.currentIndexChanged.connect(self.indexRepresChanged)

        self.SizePart = 0.3
        self.SizePartWidget = DoubleSpinBox('Size of particules', 'Plot size of particles', self.SizePart)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.SizePartWidget)

        self.NbBinsX = 300
        self.NbBinsXWidget = SpinBox('X Bining', 'Number of bins in the x-axis', self.NbBinsX, 1)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsXWidget)
        self.NbBinsXWidget.setEnabled(False)

        self.NbBinsY = 300
        self.NbBinsYWidget = SpinBox('Y Bining', 'Number of bins in the y-axis', self.NbBinsY, 1)
        self.NbBinsXWidget.Layout.addWidget(self.NbBinsYWidget)
        self.NbBinsYWidget.setEnabled(False)

        # Checkbox pour la convolution
        self.CheckConvolution = CheckBox('Simulate PSF',  'Convolve with a Gaussian kernel to simulate an instrumental PSF')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckConvolution)
        self.CheckConvolution.Layout.addSpacing(20)
        self.FWHMConv = 1
        self.FWHMConvWidget = DoubleSpinBox('Resolution', 'FWHM of the instrument [AU]', self.FWHMConv, 0, None, 1, 1)
        # self.SigmaConvWidget.Layout.insertSpacing(0, 60)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.SigmaConvWidget)
        self.CheckConvolution.Layout.addWidget(self.FWHMConvWidget)
        self.FWHMConvWidget.setEnabled(self.CheckConvolution.CheckParam.isChecked())
        self.CheckConvolution.CheckParam.stateChanged.connect(lambda: self.FWHMConvWidget.setEnabled(self.CheckConvolution.CheckParam.isChecked()))
        self.CheckConvolution.setEnabled(False)

        self.indexViewChanged(self.indexView)

    def indexViewChanged(self, value):
        self.indexView = value

        # Prevent switching to another view while in density representation
        if self.indexRepres == 1:  # If density representation is selected
            self.RepresWidget.ComboParam.setCurrentIndex(0)  # Reset to scatter representation
            self.indexRepresChanged(0)  # Update parameters accordingly

        if self.indexView == 0:
            self.RepresWidget.setEnabled(True)
            self.WidgetPlotXY.setVisible(True)
            self.WidgetPlotXZ.setVisible(False)
            self.WidgetPlotXYZ.setVisible(False)
        elif self.indexView == 1:
            self.RepresWidget.setEnabled(True)
            self.WidgetPlotXY.setVisible(False)
            self.WidgetPlotXZ.setVisible(True)
            self.WidgetPlotXYZ.setVisible(False)
        elif self.indexView == 2:
            self.RepresWidget.ComboParam.setCurrentIndex(0)
            self.RepresWidget.setEnabled(False)
            self.WidgetPlotXY.setVisible(False)
            self.WidgetPlotXZ.setVisible(False)
            self.WidgetPlotXYZ.setVisible(True)
        self.refresh_plots()

    def indexRepresChanged(self, value):
        self.indexRepres = value
        if self.indexRepres == 0:
            self.NbBinsXWidget.setEnabled(False)
            self.NbBinsYWidget.setEnabled(False)
            self.CheckConvolution.setEnabled(False)
            self.SizePartWidget.setEnabled(True)
        elif self.indexRepres == 1:
            self.NbBinsXWidget.setEnabled(True)
            self.NbBinsYWidget.setEnabled(True)
            self.CheckConvolution.setEnabled(True)
            self.SizePartWidget.setEnabled(False)
        self.refresh_plots()

    def UpdateParams(self):
        # self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        self.SizePart = self.SizePartWidget.SpinParam.value()
        self.NbBinsX = self.NbBinsXWidget.SpinParam.value()
        self.NbBinsY = self.NbBinsYWidget.SpinParam.value()
        self.FWHMConv = self.FWHMConvWidget.SpinParam.value()

    def Ellipse2(self, a, e, Ex, Ey, Ez, Epx, Epy, Epz):
        E = np.linspace(-pi, pi, 100)
        try:    
            x = Ex*a*(np.cos(E)-e) + Epx*a*sqrt(1-e**2)*np.sin(E)
            y = Ey*a*(np.cos(E)-e) + Epy*a*sqrt(1-e**2)*np.sin(E)
            z = Ez*a*(np.cos(E)-e) + Epz*a*sqrt(1-e**2)*np.sin(E)
            return x, y, z
        except Exception as e:
            print(f'Error in the ellipse computation: e > 1')
            return np.array([]), np.array([]), np.array([])

    def general_plot(self):
        # Bodies positions
        self.Xbf, self.Ybf, self.Zbf = [], [], []
        self.Xbm, self.Ybm, self.Zbm = [], [], []
        for k in range(self.NbBodies_m[self.IndexSnap]):
            # xbf, ybf, zbf = self.Ellipse(af[k][isave], ef[k][isave], i[k][isave]*pi/180, W[k][isave]*pi/180, w[k][isave]*pi/180)
            # print(self.a_m[self.IndexSnap][k], self.e_m[self.IndexSnap][k])
            xbm, ybm, zbm = self.Ellipse2(self.a_m[self.IndexSnap][k], self.e_m[self.IndexSnap][k], self.Ex[self.IndexSnap][k], self.Ey[self.IndexSnap][k], self.Ez[self.IndexSnap][k], self.Epx[self.IndexSnap][k], self.Epy[self.IndexSnap][k], self.Epz[self.IndexSnap][k])
            # self.Xbf.append(xbf)
            # self.Ybf.append(ybf)
            # self.Zbf.append(zbf)
            self.Xbm.append(xbm)
            self.Ybm.append(ybm)
            self.Zbm.append(zbm)
        
        # Specific variables
        self.NbBodies = self.NbBodies_m[self.IndexSnap]
        self.t = self.t_m[self.IndexSnap]/10**6

        # X limits
        self.LimDefault = 1.1*np.max(self.a_m[0])*(1+np.max(self.e_m[0]))

    def PlotXY(self):

        # Update of parameters
        self.try_UpdateParams(self.WidgetPlotXY)

        # General plot
        self.general_plot()

        # add subplot
        self.SubplotXY = self.WidgetPlotXY.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')
        # if self.CheckEqualAspect.CheckParam.isChecked(): self.SubplotXY.set_box_aspect(1)

        # X, Y current limits
        xlim_init = (-self.LimDefault, self.LimDefault)
        ylim_init = (-self.LimDefault, self.LimDefault)
        (Xmin, Xmax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['xlim'] if len(self.WidgetPlotXY.history)!=0 else xlim_init
        (Ymin, Ymax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['ylim'] if len(self.WidgetPlotXY.history)!=0 else ylim_init

        # Bodies' positions and orbits
        for k in range(self.NbBodies):
            if self.CheckOrbits.CheckParam.isChecked():
                self.SubplotXY.plot(self.Xbm[k], self.Ybm[k], color=self.colorList[k], linestyle='--', label="Orbit of "+str(k+1))
            if self.CheckBodies.CheckParam.isChecked():
                self.SubplotXY.plot(self.X[self.IndexSnap][k], self.Y[self.IndexSnap][k], marker='.', markersize=self.SizeBodies, color=self.colorList[k], label="Marker of "+str(k+1))

        # Plot
        if self.indexRepres == 0:
            self.SubplotXY.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Y[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0, label='Particles')
        
        elif self.indexRepres == 1:
            if self.CheckConvolution.CheckParam.isChecked():
                # Convolve positions
                x, y = self.convolve_positions(
                    x=self.X[self.IndexSnap][self.NbBodies-1:], 
                    y=self.Y[self.IndexSnap][self.NbBodies-1:], 
                    fwhm_au=self.FWHMConv,
                    extent=self.LimDefault)
            else:
                x = self.X[self.IndexSnap][self.NbBodies-1:]
                y = self.Y[self.IndexSnap][self.NbBodies-1:]
                
            hist, xedges, yedges = np.histogram2d(
                x, 
                y, 
                range=[[Xmin, Xmax], [Ymin, Ymax]], 
                bins=[self.NbBinsX, self.NbBinsY])
            norm_hist = hist
            im = self.SubplotXY.imshow(
                norm_hist.T, 
                interpolation='bicubic', 
                extent=[Xmin, Xmax, Ymin, Ymax], 
                cmap='magma', 
                origin='lower',
                label='Colormap')
            
            # Create a colorbar axis with the same height as the main plot
            divider = make_axes_locatable(self.SubplotXY)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = self.WidgetPlotXY.Canvas.fig.colorbar(im, cax=cax, ticks=[], label='Density')

        # Time
        self.SubplotXY.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.SubplotXY.transAxes, label='Time')

        # Plot features
        # self.SubplotXY.set_title('t='+str(round(self.t, 1))+' Myr')
        self.SubplotXY.set_xlabel('X [AU]')
        self.SubplotXY.set_ylabel('Y [AU]')
        self.SubplotXY.set_xlim(xlim_init)
        self.SubplotXY.set_ylim(ylim_init)


    def PlotXZ(self):

        # Update of parameters
        self.try_UpdateParams(self.WidgetPlotXZ)

        # General plot
        self.general_plot()

        # add subplot
        self.SubplotXZ = self.WidgetPlotXZ.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')
        # if self.CheckEqualAspect.CheckParam.isChecked(): self.SubplotXZ.set_box_aspect(1)

        # X, Y current limits
        xlim_init = (-self.LimDefault, self.LimDefault)
        zlim_init = (-self.LimDefault, self.LimDefault)
        (Xmin, Xmax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['xlim'] if len(self.WidgetPlotXY.history)!=0 else xlim_init
        (Zmin, Zmax) = self.WidgetPlotXY.history[self.WidgetPlotXY.history_index]['ylim'] if len(self.WidgetPlotXY.history)!=0 else zlim_init

        # Plot
        for k in range(self.NbBodies):
            if self.CheckOrbits.CheckParam.isChecked():
                self.SubplotXZ.plot(self.Xbm[k], self.Zbm[k], color=self.colorList[k], linestyle='--', label="Orbit of "+str(k+1))
            if self.CheckBodies.CheckParam.isChecked():
                self.SubplotXZ.plot(self.X[self.IndexSnap][k], self.Z[self.IndexSnap][k], marker='.', markersize=self.SizeBodies, color=self.colorList[k], label="Marker of "+str(k+1))

        if self.indexRepres == 0:
            self.SubplotXZ.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Z[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0, label='Particles')
        
        elif self.indexRepres == 1:
            if self.CheckConvolution.CheckParam.isChecked():
                # Convolve positions
                x, z = self.convolve_positions(
                    x=self.X[self.IndexSnap][self.NbBodies-1:], 
                    y=self.Z[self.IndexSnap][self.NbBodies-1:], 
                    fwhm_au=self.FWHMConv,
                    extent=self.LimDefault)
            else:
                x = self.X[self.IndexSnap][self.NbBodies-1:]
                z = self.Z[self.IndexSnap][self.NbBodies-1:]
            
            hist, xedges, yedges = np.histogram2d(
                x, 
                z, 
                range=[[Xmin, Xmax], [Zmin, Zmax]], 
                bins=[self.NbBinsX, self.NbBinsY])
            norm_hist = hist
            im = self.SubplotXZ.imshow(
                norm_hist.T, 
                interpolation='bicubic', 
                extent=[Xmin, Xmax, Zmin, Zmax], 
                cmap='magma', 
                origin='lower',
                label='Colormap')
            
            # Create a colorbar axis with the same height as the main plot
            divider = make_axes_locatable(self.SubplotXZ)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = self.WidgetPlotXZ.Canvas.fig.colorbar(im, cax=cax, ticks=[], label='Density')
       
        # Time
        self.SubplotXZ.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.SubplotXZ.transAxes, label='Time')

        # Plot features
        self.SubplotXZ.set_xlabel('X [AU]')
        self.SubplotXZ.set_ylabel('Z [AU]')
        self.SubplotXZ.set_xlim(xlim_init)
        self.SubplotXZ.set_ylim(zlim_init)
        

    def PlotXYZ(self):

        # Update of parameters
        self.try_UpdateParams(self.WidgetPlotXYZ)

        # General plot
        self.general_plot()

        # Add subplot
        self.SubplotXYZ = self.WidgetPlotXYZ.Canvas.fig.add_subplot(111, aspect='equal', projection='3d', label='Main plot')
        # if self.CheckEqualAspect.CheckParam.isChecked(): self.SubplotXYZ.set_box_aspect(1)

        # X, Y, Z limits
        xlim_init = [-self.LimDefault, self.LimDefault]
        ylim_init = [-self.LimDefault, self.LimDefault]
        zlim_init = [-self.LimDefault, self.LimDefault]

        # Plot
        for k in range(self.NbBodies):
            if self.CheckOrbits.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.Xbm[k], self.Ybm[k], self.Zbm[k], color=self.colorList[k], linestyle='--', label=f"Orbit of {k+1}")
            if self.CheckBodies.CheckParam.isChecked():
                self.SubplotXYZ.plot(self.X[self.IndexSnap][k], self.Y[self.IndexSnap][k], self.Z[self.IndexSnap][k], 
                                     marker='.', markersize=self.SizeBodies, color=self.colorList[k], label=f"Marker of {k+1}")

        self.SubplotXYZ.scatter(self.X[self.IndexSnap][self.NbBodies-1:], 
                                self.Y[self.IndexSnap][self.NbBodies-1:], 
                                self.Z[self.IndexSnap][self.NbBodies-1:], 
                                s=self.SizePart, c='black', linewidth=0, label='Particles')
        
        # Time
        self.SubplotXYZ.text2D(1.01, 1.01, 't='+str(round(self.t, 1))+' Myr', 
                               horizontalalignment='right', verticalalignment='top', 
                               transform=self.SubplotXYZ.transAxes, label='Time')

        # Ensure axis limits are finite
        self.SubplotXYZ.set_xlim(xlim_init)
        self.SubplotXYZ.set_ylim(ylim_init)
        self.SubplotXYZ.set_zlim(zlim_init)

        # Plot features
        self.SubplotXYZ.set_xlabel('X [AU]')
        self.SubplotXYZ.set_ylabel('Y [AU]')
        self.SubplotXYZ.set_zlabel('Z [AU]')


class RadProfile(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m, X, Y, Z):
        super().__init__('Radial profile', "Particules' integrated radial profile")

        # Data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.e_m = e_m
        self.X = X
        self.Y = Y
        self.Z = Z 

        # Beta Pic other curves
        # self.profileAug = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Augereau_profile.dat', dtype = float))
        # self.profileNE = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Dent_profile_NE.dat', dtype = float))
        # self.profileSW = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Dent_profile_SW.dat', dtype = float))

        # Window plot initialisation
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylabel=False)

        # Parameters initialisation
        self.LimDefault = 1.1*np.max(self.a_m[0])*(1+np.max(self.e_m[0]))
        self.InitParams()

    def InitParams(self):
        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

        # Histogram
        self.Ordinate = ComboBox('Ordinate', 'Choice of ordinate', ['Number of particules', 'Surface density', 'Other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.Ordinate)
        # self.Ordinate.ComboParam.currentIndexChanged.connect(self.WidgetPlot.reset_plots)
        self.Ordinate.ComboParam.currentIndexChanged.connect(self.WidgetPlot.refresh_plot)

        self.FormulaTextEdit = LineEdit(None, 'Custom formula for the ordinate (use "r" and "counts" for radius and number of particles in each bin)', None)
        self.FormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.FormulaTextEdit.setVisible(False)
        self.Ordinate.ComboParam.currentIndexChanged.connect(lambda: self.FormulaTextEdit.setVisible(self.Ordinate.ComboParam.currentText() == 'Other'))
        self.WindowPlot.WidgetParam.Layout.addWidget(self.FormulaTextEdit)
        # self.FormulaTextEdit.EditParam.textChanged.connect(self.WidgetPlot.refresh_plot)

        self.Norm = ComboBox('Normalisation', 'Choice of the normalisation', ['None', 'max=1', 'sum=1'])
        self.Norm.ComboParam.setCurrentIndex(0)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.Norm)
        # self.Norm.ComboParam.currentIndexChanged.connect(self.WidgetPlot.reset_plots)
        self.Norm.ComboParam.currentIndexChanged.connect(self.WidgetPlot.refresh_plot)

        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Bining', 'Number of bins', self.NbBins, 1)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)
        # self.NbBinsWidget.Layout.addSpacing(20)
        # self.BinSize = 0.1
        # self.BinSizeWidget = DoubleSpinBox('Bin size', 'Size of the bins [AU]', self.BinSize, 0.01, None, 0.01, 0.01)
        # self.NbBinsWidget.Layout.addWidget(self.BinSizeWidget)

        # Checkbox pour l'interpolation
        self.CheckInterpolation = CheckBox('Interpolate', 'Interpolate the radial profile')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckInterpolation)

        # Checkbox pour la convolution
        self.CheckConvolution = CheckBox('Simulate PSF',  'Convolve with a gaussian kernel to simulate an instrumental PSF')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckConvolution)
        self.CheckConvolution.Layout.addSpacing(20)
        self.FWHMConv = 1
        self.FWHMConvWidget = DoubleSpinBox('Resolution', 'FWHM of the instrument [AU]', self.FWHMConv, 0, 100, 1, 1)
        # self.SigmaConvWidget.Layout.insertSpacing(0, 60)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.SigmaConvWidget)
        self.CheckConvolution.Layout.addWidget(self.FWHMConvWidget)
        self.FWHMConvWidget.setEnabled(self.CheckConvolution.CheckParam.isChecked())
        self.CheckConvolution.CheckParam.stateChanged.connect(lambda: self.FWHMConvWidget.setEnabled(self.CheckConvolution.CheckParam.isChecked()))

        # Checkbox pour la detection de bord
        self.CheckEdge = CheckBox('Edge detection', 'Detect edge in the radial profile')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEdge)
        self.CheckEdge.Layout.addSpacing(20)
        self.LowerBound = 0
        self.LowerBoundWidget = DoubleSpinBox('Prior', 'Lower bound of the edge [AU]', self.LowerBound, 0, self.LimDefault, 1, 1)
        self.CheckEdge.Layout.addWidget(self.LowerBoundWidget)
        self.LowerBoundWidget.setEnabled(self.CheckEdge.CheckParam.isChecked())
        self.CheckEdge.CheckParam.stateChanged.connect(lambda: self.LowerBoundWidget.setEnabled(self.CheckEdge.CheckParam.isChecked()))
        # self.LabelInf = QLabel('<')
        # self.CheckEdge.Layout.addWidget(self.LabelInf)
        # self.LabelInf.setEnabled(self.CheckEdge.CheckParam.isChecked())
        # self.CheckEdge.CheckParam.stateChanged.connect(lambda: self.LabelInf.setEnabled(self.CheckEdge.CheckParam.isChecked()))
        # self.Edge = 0.00
        # self.EdgeWidget = DoubleSpinBox(None, 'Prior edge [AU]', self.Edge, 0, self.LimDefault, 1, 2)
        # self.CheckEdge.Layout.addWidget(self.EdgeWidget)
        # self.EdgeWidget.setEnabled(self.CheckEdge.CheckParam.isChecked())
        # self.CheckEdge.CheckParam.stateChanged.connect(lambda: self.EdgeWidget.setEnabled(self.CheckEdge.CheckParam.isChecked()))
        self.LabelInterval = QLabel('<->')
        self.CheckEdge.Layout.addWidget(self.LabelInterval)
        self.LabelInterval.setEnabled(self.CheckEdge.CheckParam.isChecked())
        self.CheckEdge.CheckParam.stateChanged.connect(lambda: self.LabelInterval.setEnabled(self.CheckEdge.CheckParam.isChecked()))
        self.UpperBound = self.LimDefault
        self.UpperBoundWidget = DoubleSpinBox(None, 'Upper bound of the edge [AU]', self.UpperBound, 0, self.LimDefault, 1, 1)
        self.CheckEdge.Layout.addWidget(self.UpperBoundWidget)
        self.UpperBoundWidget.setEnabled(self.CheckEdge.CheckParam.isChecked())
        self.CheckEdge.CheckParam.stateChanged.connect(lambda: self.UpperBoundWidget.setEnabled(self.CheckEdge.CheckParam.isChecked()))

        # Bodies' positions
        self.CheckBodies = CheckBox('Bodies position')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBodies)

        self.SizeBodies = 15
        # self.SizeBodiesWidget = SpinBox('Size of bodies', 'Plot size of bodies', self.SizeBodies)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.SizeBodiesWidget)
        
        # self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked())
        # self.CheckBodies.CheckParam.stateChanged.connect(lambda: self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked()))

        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Other curves
        self.CheckCurves = QCheckBox('Other curves')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckCurves, alignment=Qt.AlignmentFlag.AlignCenter)
        self.CheckCurves.stateChanged.connect(self.CurveStateChange)

        # self.CheckAugWidget = QCheckBox('Aug 2001')
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckAugWidget)
        # self.CheckAugWidget.setEnabled(self.CheckCurves.isChecked())

        # self.CheckDentWidget = QCheckBox('Dent 2014')
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckDentWidget)
        # self.CheckDentWidget.setEnabled(self.CheckCurves.isChecked())

        self.ButAddCurve = QPushButton('+')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ButAddCurve)
        self.ButAddCurve.clicked.connect(self.AddCurveWidget)
        self.ButAddCurve.setEnabled(self.CheckCurves.isChecked())

        self.CurveWidgets = []
        self.CurveWidgetsId = []
        self.NbWidgets0 = self.WindowPlot.WidgetParam.Layout.count()
        self.c = 0 # counter



    def AddCurveWidget(self):
        self.CurveWidget = CurveClass()
        self.CurveWidgets.append(self.CurveWidget)
        self.c += 1
        self.CurveWidget.Id = self.c
        self.CurveWidgetsId.append(self.CurveWidget.Id)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CurveWidget)
        self.CurveWidget.SignalDel.connect(self.DelCurveWidget)

    def DelCurveWidget(self, Id):
        index = self.CurveWidgetsId.index(Id)
        self.CurveWidgets.pop(index)
        self.CurveWidgetsId.pop(index)

        for i in reversed(range(self.NbWidgets0, self.WindowPlot.WidgetParam.Layout.count())): 
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        
        for i in range(len(self.CurveWidgets)):
            self.WindowPlot.WidgetParam.Layout.addWidget(self.CurveWidgets[i])

        self.WindowPlot.WidgetParam.resize(self.WindowPlot.WidgetParam.minimumSize())

    def CurveStateChange(self):
        CurveState = self.CheckCurves.isChecked()
        # self.CheckAugWidget.setEnabled(CurveState)
        # self.CheckDentWidget.setEnabled(CurveState)
        self.ButAddCurve.setEnabled(CurveState)
        for i in range(len(self.CurveWidgets)):
            self.CurveWidgets[i].setEnabled(CurveState)
    
    def UpdateParams(self):
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        self.FWHMConv = self.FWHMConvWidget.SpinParam.value()
        self.LowerBound = self.LowerBoundWidget.SpinParam.value()
        # self.Edge = self.EdgeWidget.SpinParam.value()
        self.UpperBound = self.UpperBoundWidget.SpinParam.value()
        # self.EdgeSigma = float(self.EdgeSigmaWidget.text().split(' ')[-1])
        # self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        if self.CheckCurves.isChecked():
            self.CurvePaths = []
            self.CurveLabels = []
            for i in range(len(self.CurveWidgets)):
                self.CurvePaths.append(self.CurveWidgets[i].PathFile)
                name = self.CurveWidgets[i].NameFile.EditParam.text()
                name_without_ext = os.path.splitext(name)[0]
                self.CurveLabels.append(name_without_ext)

    # def Rebin(self, data, resolution):    
    #     """
    #     Rebin the data to a specified resolution.
    #     """
    #     smooth_data, bins_edges = histogram(data, bins=(np.max(data)-np.min(data))/resolution)
    #     return smooth_data 

    # Plot
    def Plot(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

        # Radius computation
        if self.CheckConvolution.CheckParam.isChecked():
            # Convolve positions
            x_conv, y_conv = self.convolve_positions(
                x=self.X[self.IndexSnap], 
                y=self.Y[self.IndexSnap], 
                fwhm_au=self.FWHMConv,
                extent=self.LimDefault)
            self.R = np.sqrt(x_conv**2 + y_conv**2)
        else:
            self.R = np.sqrt(self.X[self.IndexSnap]**2+self.Y[self.IndexSnap]**2)

        # X, Z limits
        xlim_init = (0, self.LimDefault)
        # xlim_init = (0, self.LimDefault)
        xlim = self.WidgetPlot.history[self.WidgetPlot.history_index]['xlim'] if self.WidgetPlot.history else xlim_init
        Rmin, Rmax = xlim[0], xlim[1]
        
        # Specific variables
        self.t = self.t_m[self.IndexSnap]/10**6

        # Plot with current parameters
        if self.CheckBodies.CheckParam.isChecked():
            for k in range(self.NbBodies_m[self.IndexSnap]):
                self.Subplot.plot(self.a_m[self.IndexSnap][k], 0, color=self.colorList[k], marker='.', markersize=self.SizeBodies)
    
        histCount, histX = np.histogram([x for x in self.R if Rmin < x < Rmax], bins=self.NbBins)

        # Surface density computation
        r = (histX[:-1] + histX[1:]) / 2  # midpoints of bins
        if self.Ordinate.ComboParam.currentIndex() == 1:
            dr = np.diff(histX)  # bin widths
            area = 2 * pi * r * dr
            # Avoid division by zero for r=0
            area[area == 0] = np.nan
            histCount = histCount / area  # true surface density
        elif self.Ordinate.ComboParam.currentText() == 'Other':
            counts = histCount
            formula = self.FormulaTextEdit.EditParam.text()
            try:
                histCount = eval(formula)
                np.shape(histCount) == np.shape(counts)  # Check if histCount is valid
            except Exception as e:
                print(f"Error in custom formula: {e}")
                histCount = counts
                
        # Normalisation computation
        self.NormDiv = 1
        if self.Norm.ComboParam.currentIndex()==1: self.NormDiv = np.max(histCount) # normalisation of max equal one
        elif self.Norm.ComboParam.currentIndex()==2: self.NormDiv = np.sum(histCount) # normalisation of sum equal one
        histCount = histCount/self.NormDiv

        # Stairs
        self.Subplot.stairs(histCount, histX, label='Simulation', linewidth=1, color='black')

        # Interpolation si la checkbox est cochée
        if self.CheckInterpolation.CheckParam.isChecked():
            try:
                # Utilisation des milieux des bins pour l'interpolation
                histX_mid = (histX[:-1] + histX[1:]) / 2  # Calcul des milieux des bins
                interp_func = interp1d(histX_mid, histCount, kind='cubic', fill_value="extrapolate")
                
                # Génération de points interpolés
                histX_interp = np.linspace(histX_mid[0], histX_mid[-1], 1000)  # Plus de points pour une courbe lisse
                histCount_interp = interp_func(histX_interp)
                
                # Tracé de la courbe interpolée
                self.Subplot.plot(histX_interp, histCount_interp, label='Interpolation', linewidth=1, color='green')
            except Exception as e:
                print(f"Interpolation failed: {e}")

        # Fenêtre autour du bord prior pour le fit
        if self.LowerBound > self.UpperBound:
            self.LowerBound, self.UpperBound = self.UpperBound, self.LowerBound
        if self.LowerBound < self.Subplot.get_xlim()[0]:
            self.LowerBound = self.Subplot.get_xlim()[0]
        if self.UpperBound > self.Subplot.get_xlim()[1]:
            self.UpperBound = self.Subplot.get_xlim()[1]

        self.LowerBoundWidget.SpinParam.setValue(self.LowerBound)
        self.UpperBoundWidget.SpinParam.setValue(self.UpperBound)

        # Edge detection si la checkbox est cochée
        if self.CheckEdge.CheckParam.isChecked():
            try:
                # Calcul de la dérivée (gradient)
                smooth_hist = gaussian_filter1d(histCount, sigma=2)  # Lissage avant dérivation
                gradient = np.gradient(smooth_hist, histX[:-1])

                mask = (self.LowerBound < histX[:-1]) & (histX[:-1] < self.UpperBound)
                histX_mask = histX[:-1][mask]
                gradient_mask = gradient[mask]

                # Détection du bord prior
                # if self.Edge !=0 and self.LowerBound < self.Edge < self.UpperBound:
                #     # print(f"Using user-defined edge: {self.Edge}")
                #     edge_prior_index = np.argmin(np.abs(histX[:-1] - self.Edge))  # Trouver l'index le plus proche de l'Edge
                #     edge_prior = self.Edge  # Utilisation de la valeur de bord définie par l'utilisateur
                # else:
                    # print("No user-defined edge, finding edge from gradient")
                edge_prior_index = np.argmax(np.abs(gradient_mask))
                edge_prior = histX_mask[edge_prior_index]

                gradient_sign = np.sign(gradient_mask[edge_prior_index])

                # # Fenêtre autour du bord prior pour le fit
                # fit_mask = (self.LowerBound < histX[:-1]) & (histX[:-1] < self.UpperBound)
                # x_fit = histX[:-1][fit_mask]
                # y_fit = gradient[fit_mask]

                # Fit d'une gaussienne
                if gradient_sign > 0:
                    p0 = [np.max(gradient_mask), edge_prior, 5]  # amplitude, centre, largeur initiale
                    bounds = ([0, self.LowerBound, 0], [np.inf, self.UpperBound, np.inf])  # Ajuster les bornes si nécessaire
                else:
                    p0 = [np.min(gradient_mask), edge_prior, 5]
                    bounds = ([-np.inf, self.LowerBound, 0], [0, self.UpperBound, np.inf])  # Ajuster les bornes si nécessaire
                popt, _ = curve_fit(self.gaussian, histX_mask, gradient_mask, p0=p0, bounds=bounds)

                # Calcul de la position du bord et de l'écart-type
                edge_position = popt[1]
                edge_sigma = popt[2]
                # print(f"Edge detected at {edge_position} with sigma {edge_sigma}")
                lower_sigma = edge_position - edge_sigma
                upper_sigma = edge_position + edge_sigma

                # Tracé du point de bord
                self.Subplot.axvline(edge_position, color='red', label='Edge detected', linewidth=1)
                self.Subplot.text(edge_position, 0.83 * self.Subplot.get_ylim()[1], s='{}'.format(np.around(edge_position, 2)), color='red', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red'), fontsize=9, horizontalalignment='center', verticalalignment='bottom')

                self.Subplot.axvline(edge_position + edge_sigma, color='red', linestyle='--')
                self.Subplot.text(lower_sigma, 0.8 * self.Subplot.get_ylim()[1], s='-{}'.format(np.around(edge_sigma, 2)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red'), fontsize=8, horizontalalignment='right', verticalalignment='top')

                self.Subplot.axvline(edge_position - edge_sigma, color='red', linestyle='--')
                self.Subplot.text(upper_sigma, 0.8 * self.Subplot.get_ylim()[1], s='+{}'.format(np.around(edge_sigma, 2)), color='r', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red'), fontsize=8, horizontalalignment='left', verticalalignment='top')      

            except Exception as e:
                print(f"Edge detection failed: {e}")

        # Other curves
        # if self.CheckAugWidget.isChecked(): self.Subplot.plot(self.profileAug[0], self.profileAug[1], color='blue', linestyle='dashed', linewidth=0.5, label='Aug+2001')
        # if self.CheckDentWidget.isChecked(): self.Subplot.plot((self.profileNE[0]+self.profileSW[0])/2, (self.profileNE[1]+self.profileSW[1])/2, color='red', linestyle='dashed', linewidth=0.5, label='Dent+2014')
        
        if self.CheckCurves.isChecked(): 
            try:
                for i in range(len(self.CurveWidgets)):
                    # Automatic delimiter detection
                    delimiter = find_delimiter(self.CurvePaths[i])
                    # print(delimiter)
                    self.Curve = np.transpose(np.genfromtxt(self.CurvePaths[i], delimiter=delimiter, dtype=float))
                    self.Subplot.plot(self.Curve[0], self.Curve[1], linestyle='dashed', linewidth=1, label=self.CurveLabels[i])
                    self.Subplot.legend()
            except Exception as e:
                print(f'Wrong format: {e}')
                print('Two columns of equal dimensions are required, corresponding to the x and y points of the curve')

        # Time
        self.Subplot.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.Subplot.transAxes)

        # Plot features
        # if self.CheckCurves.isChecked() and len(self.CurvePaths)>0 : self.Subplot.legend()
        self.Subplot.set_xlabel('Radius [AU]')
        if self.Ordinate.ComboParam.currentIndex()==0: self.Subplot.set_ylabel('Number of particules')
        elif self.Ordinate.ComboParam.currentIndex()==1: self.Subplot.set_ylabel('Surface density [arbitrary unit]')
        elif self.Ordinate.ComboParam.currentText() == 'Other': self.Subplot.set_ylabel('Custom ordinate [arbitrary unit]')

    def gaussian(self, x, a, mu, sigma):
        return a * np.exp(-(x - mu)**2 / (2 * sigma**2))


class DiagramAE(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m):
        super().__init__('Diagram a=f(e)', "Bodies' evolution of semi-major axis as a function of excentricity")

        # General data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.e_m = e_m

        # Plot initialisation
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

        # Parameters initialisation
        self.InitParams()

    def InitParams(self):
        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

        # Bodies' positions
        self.CheckBodies = CheckBox("Bodies position")
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBodies)
        self.CheckBodies.CheckParam.checkStateChanged.connect(self.WidgetPlot.refresh_plot)

        self.SizeBodies = 15
        # self.SizeBodiesWidget = SpinBox('Size of bodies', 'Plot size of bodies', self.SizeBodies)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.SizeBodiesWidget)
        
        # self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked())
        # self.CheckBodies.CheckParam.stateChanged.connect(lambda: self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked()))

        # self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Particules' size
        self.SizePart = 0.05
        self.SizePartWidget = DoubleSpinBox('Size of particules', 'Plot size of particles', self.SizePart, ParamPrecision=2, ParamIncrement=0.01)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.SizePartWidget)

        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Orbital resonance
        self.CheckRes = QCheckBox('Orbital resonance')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckRes, alignment=Qt.AlignmentFlag.AlignCenter)
        self.CheckRes.stateChanged.connect(self.ResStateChange)

        self.ButAddRes = QPushButton('+')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ButAddRes)
        self.ButAddRes.clicked.connect(self.AddResWidget)
        self.ButAddRes.setEnabled(False)

        self.ResWidgets = []
        self.ResWidgetsId = []
        self.NbWidgets0 = self.WindowPlot.WidgetParam.Layout.count()

        self.c = 0 # counter

    def AddResWidget(self):
        self.ResWidget = ResClass(self.NbBodies_m[self.IndexSnap])
        self.ResWidgets.append(self.ResWidget)
        self.c += 1
        self.ResWidget.Id = self.c
        self.ResWidgetsId.append(self.ResWidget.Id)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ResWidget)
        self.ResWidget.SignalDel.connect(self.DelResWidget)

    def DelResWidget(self, Id):
        index = self.ResWidgetsId.index(Id)
        self.ResWidgets.pop(index)
        self.ResWidgetsId.pop(index)

        for i in reversed(range(self.NbWidgets0, self.WindowPlot.WidgetParam.Layout.count())): 
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        
        for i in range(len(self.ResWidgets)):
            self.WindowPlot.WidgetParam.Layout.addWidget(self.ResWidgets[i])

        self.WindowPlot.WidgetParam.resize(self.WindowPlot.WidgetParam.minimumSize())

    def ResStateChange(self):
        ResState = self.CheckRes.isChecked()
        self.ButAddRes.setEnabled(ResState)
        for i in range(len(self.ResWidgets)):
            self.ResWidgets[i].setEnabled(ResState)

    # Update of parameters
    def UpdateParams(self):
        self.SizePart = self.SizePartWidget.SpinParam.value()
        # self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        if self.CheckRes.isChecked():
            self.nRefValues = []
            self.PResValues = []
            self.PRefValues = []
            for i in range(len(self.ResWidgets)):
                self.nRefValues.append(int(self.ResWidgets[i].nRefWidget.ComboParam.currentText()))
                self.PResValues.append(int(self.ResWidgets[i].PResWidget.SpinParam.value()))
                self.PRefValues.append(int(self.ResWidgets[i].PRefWidget.SpinParam.value()))
        
    # Plot
    def Plot(self):

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)
        
        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)
    
        # Specific variables
        self.a = self.a_m[self.IndexSnap][:]
        self.e = self.e_m[self.IndexSnap][:]
        self.t = self.t_m[self.IndexSnap]/10**6
        if self.CheckRes.isChecked():
            self.aRefValues = []
            self.aResValues = []
            for i in range(len(self.ResWidgets)):
                self.aRefValues.append(self.a[self.nRefValues[i]])
                self.aResValues.append((self.PResValues[i]/self.PRefValues[i])**(2/3)*self.aRefValues[i])

        # X, Y limits
        xlim_init = (0, np.max(self.a)+0.1*np.max(self.a))
        ylim_init = (0, 1)
        (Ymin, Ymax) = self.WidgetPlot.history[self.WidgetPlot.history_index]['ylim'] if self.WidgetPlot.history else ylim_init
        
        # Plot with current parameters
        if self.CheckBodies.CheckParam.isChecked():
            for k in range(self.NbBodies_m[self.IndexSnap]):
                self.Subplot.plot(self.a[k], self.e[k], color=self.colorList[k], marker='.', markersize=self.SizeBodies, label="Marker of "+str(k+1))
        self.Subplot.scatter(self.a[self.NbBodies_m[self.IndexSnap]:], self.e[self.NbBodies_m[self.IndexSnap]:], s=self.SizePart, c='black', linewidths=0)

        if self.CheckRes.isChecked():
            for i in range(len(self.ResWidgets)):
                colorRef = self.colorList[self.nRefValues[i]]
                self.Subplot.axvline(self.aResValues[i], linewidth=1, color=colorRef, linestyle='--', label='Resonance '+str(self.PResValues[i])+':'+str(self.PRefValues[i])+' of '+str(self.nRefValues[i]))
                # self.Subplot.annotate(str(self.PResValues[i])+':'+str(self.PRefValues[i]), xy=(self.aResValues[i], 0.9), bbox=dict(boxstyle='round', facecolor='white', edgecolor=colorRef))
                self.Subplot.text(self.aResValues[i], 0.9*Ymax, ' '+str(self.PResValues[i])+':'+str(self.PRefValues[i])+' ', rotation=90, color=colorRef, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=colorRef), fontsize=7, va='center', ha='center')

        # Time
        self.Subplot.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.Subplot.transAxes)

        # Plot features
        # self.Subplot.set_title('t='+str(round(self.t, 1))+' Myr')
        self.Subplot.set_xlabel('Semi-major axis [AU]')
        self.Subplot.set_xlim(xlim_init)
        self.Subplot.set_ylabel('Eccentricity')
        self.Subplot.set_ylim(ylim_init)
            

class DiagramTY(GeneralToolClass):
    def __init__(self, NbBodies_f, t_f, a_f, e_f, i, W, w, M):
        super().__init__('Diagram y=f(t)', "Planets' evolution of an orbital parameter as a function of time")

        # General data
        self.NbBodies_f = NbBodies_f
        self.t = t_f
        self.a = a_f
        self.e = e_f
        self.i = i
        self.W = W
        self.w = w
        self.M = M

        # Plot initialisation
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylabel=False)

        # Parameters initialisation
        self.InitParams()


    def InitParams(self):
        # self.CheckTitle.setEnabled(False)

        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['a','e','i','W','w','M','irel','other'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # Orbit number
        self.ListBody = [str(k + 1) for k in range(self.NbBodies_f-1)]
        self.nOrbitWidget = ComboBox(None,'Number of the orbit counting outward', self.ListBody)
        self.ParamOrbitWidget.Layout.addWidget(self.nOrbitWidget)
        self.nOrbitWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)

        # i relatif between 2 bodies
        self.IrelWidget = QWidget()
        self.IrelLayout = QHBoxLayout()
        self.nOrbitRelWidget = ComboBox('Relative body', 'Relative body to compare with the reference', self.ListBody)
        self.IrelLayout.addWidget(self.nOrbitRelWidget)
        self.nOrbitRelWidget.ComboParam.currentIndexChanged.connect(self.reset_plots)
        self.IrelWidget.setLayout(self.IrelLayout)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.IrelWidget)
        self.IrelWidget.setVisible(False)

        # TextEdit for general formula
        self.FormulaTextEdit = LineEdit(None, 'Only variables P, a, e, i, w, W, M with [n] for orbit number and usual mathematical functions', None)
        self.FormulaTextEdit.EditParam.setPlaceholderText("Enter your formula here")
        self.FormulaTextEdit.setVisible(False)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.FormulaTextEdit)
        self.FormulaTextEdit.EditParam.textChanged.connect(self.reset_plots)

        # Connect ComboBox change event
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.FormulaTextEdit.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'other'))
        self.ParamOrbitWidget.ComboParam.currentIndexChanged.connect(lambda: self.IrelWidget.setVisible(self.ParamOrbitWidget.ComboParam.currentText() == 'irel'))

        # # Time limits
        # self.Tmin = 0
        # self.TminWidget = SpinBox('Tmin', 'Time minimum [yr]', self.Tmin)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.TminWidget)

        # self.Tmax = round(np.max(self.t_f))
        # self.TmaxWidget = SpinBox('Tmax', 'Time maximum [yr]', self.Tmax)
        # self.TminWidget.Layout.addWidget(self.TmaxWidget) 

    def UpdateParams(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.nOrbit = int(self.nOrbitWidget.ComboParam.currentText())
        self.nOrbitRel = int(self.nOrbitRelWidget.ComboParam.currentText())
        self.EvalParamOrbit = self.evaluate_ParamOrbit('self.')
        # self.Tmin = self.TminWidget.SpinParam.value()
        # self.Tmax = self.TmaxWidget.SpinParam.value()

    # def ToggleFormulaTextEdit(self):
    #     """Toggle the visibility of the formula text edit based on the ComboBox selection."""
    #     if self.ParamOrbitWidget.ComboParam.currentIndex() == self.ParamOrbitWidget.ComboParam.count() - 1:
    #         self.FormulaTextEdit.setVisible(True)
    #     else:
    #         self.FormulaTextEdit.setVisible(False)
            
    def evaluate_ParamOrbit(self, prefixe):
        """Evaluate the parameter orbit based on the current widget values."""
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        if self.ParamOrbit == 'irel' or self.ParamOrbit == 'other':
            if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
                if self.nOrbitRel == self.nOrbit:
                    print('nBodyRel can not be the same as nBody')
                    return None
                else:
                    formula = f'arccos(cos(i[{self.nOrbit}])*cos(i[{self.nOrbitRel}])+cos(W[{self.nOrbit}]-W[{self.nOrbitRel}])*sin(i[{self.nOrbit}])*sin(i[{self.nOrbitRel}]))'
            elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
                formula = self.FormulaTextEdit.EditParam.text()
            return self.evaluate_formula(formula, prefixe, self.nOrbit)
        else:
            # self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
            # print('not irel or other')
            return eval(f'{prefixe}{self.ParamOrbit}')[self.nOrbit]

    # Plot
    def Plot(self):
        
        # Subplot initialization
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

        # Update parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Check if the parameters are valid for plotting
        if self.EvalParamOrbit is None:
            print('No data to plot or variance is zero.')
            self.Subplot.figure.canvas.draw()
            return

        # try:
        #     self.UpdateParams()
        # except Exception as e:
        #     print('Wrong Parameters: ', e)
        #     self.Subplot.figure.canvas.draw()
        
        # # Specific data
        # t = [self.t_f[self.nOrbit],'Time [yr]']
        # a = [self.a_f[self.nOrbit],'Semi-major axis [AU]']
        # e = [self.e_f[self.nOrbit],'Eccentricity']
        # i = [self.i[self.nOrbit],'Inclinaison [deg]']
        # W = [self.W[self.nOrbit],'Longitude of ascending node [deg]']
        # w = [self.w[self.nOrbit],'Argument of periastron [deg]']
        # M = [self.M[self.nOrbit],'Initial mean longitude']
        
        # Plot with current parameters
        # self.EvalParamOrbit = eval(self.ParamOrbit)
        self.Subplot.plot(self.t[self.nOrbit], self.EvalParamOrbit, linewidth=0.5)

        # Plot features
        self.Subplot.set_xlabel('Time [yr]')
        # self.Subplot.set_xlim(self.Tmin, self.Tmax)
        # self.Subplot.set_ylabel(self.EvalParamOrbit[1])
        if self.ParamOrbitWidget.ComboParam.currentText() == 'irel':
            ylabel_init = r'i$_{rel}$ [°]'
        elif self.ParamOrbitWidget.ComboParam.currentText() == 'other':
            ylabel_init = self.FormulaTextEdit.EditParam.text()
        else:
            ylabel_init = self.LabelOf(self.ParamOrbit)+' '+self.UnitOf(self.ParamOrbit)
        self.Subplot.set_ylabel(ylabel_init)
        
        # # Update canvas
        # self.WindowPlot.WidgetPlots.Canvas.draw()


class DiagramXY(GeneralToolClass):
    def __init__(self, NbBodies_f, t_f, a_f, e_f, i, W, w, M):
        super().__init__('Diagram y=f(x)', "Planets' evolution of an orbital parameter Y as a function of an other X")

        # General data
        self.NbBodies_f = NbBodies_f
        self.t = t_f
        self.a = a_f
        self.e = e_f
        self.i = i
        self.W = W
        self.w = w
        self.M = M

        # Plot initialisation
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, xlim=True, ylim=True)

        # Parameters initialisation
        self.InitParams()


    def InitParams(self):
        # Equal aspect ratio
        self.CheckEqualAspect = CheckBox('Equal aspect ratio', 'Force equal aspect ratio for the plot')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckEqualAspect)
        
        # Ordinate quantity formula
        self.YFormula = ''
        self.YFormulaWidget = LineEdit('Y formula','Ordinate quantity by combining t, a, e, i, w, W, with [n] for orbit number counting outwards, and usual mathematical functions', self.YFormula)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YFormulaWidget)
        self.YFormulaWidget.setMinimumWidth(300)
        self.YFormulaWidget.EditParam.textChanged.connect(self.reset_plots)

        # Abscissa quantity formula
        self.XFormula = '' 
        self.XFormulaWidget = LineEdit('X formula','Abscissa quantity by combining t, a, e, i, w, W, with [n] for orbit number counting outwards, and usual mathematical functions', self.XFormula)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XFormulaWidget)
        self.XFormulaWidget.setMinimumWidth(300)
        self.XFormulaWidget.EditParam.textChanged.connect(self.reset_plots)



    def UpdateParams(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        self.YFormula = self.YFormulaWidget.EditParam.text()
        self.XFormula = self.XFormulaWidget.EditParam.text()
        self.Y = self.evaluate_formula(self.YFormula, 'self.')
        self.X = self.evaluate_formula(self.XFormula, 'self.')
        print(f'X: {self.X}, Y: {self.Y}')  # Debugging output
                
    # Plot
    def Plot(self):

        # Update of parameters
        self.try_UpdateParams(self.WidgetPlot)

        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)
        if self.CheckEqualAspect.CheckParam.isChecked(): self.Subplot.set_box_aspect(1)

        if self.X is None or self.Y is None or np.shape(self.X) != np.shape(self.Y):
            print('No data to plot or data shapes do not match.')
            self.Subplot.figure.canvas.draw()
            return

        self.Subplot.plot(self.X, self.Y, linewidth=0.5, label='Y=f(X)')

        # Plot features
        self.Subplot.set_xlabel(self.XFormula)
        self.Subplot.set_ylabel(self.YFormula)






### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    n=1
    PathFollowbodies = f'/Users/lacquema/Documents/Swiftdata/BetaPic/simu_bpicbcd_{n}/followbodies.dat'
    PathMextract = f'/Users/lacquema/Documents/Swiftdata/BetaPic/simu_bpicbcd_{n}/mextract.dat'       
    NbSteps, NbBodies_f, t_f, a_f, e_f, i, W, w, M = TransferDataClass.OpenFollowbodies(PathFollowbodies)
    NbSnapshots, t_m, NbBodies_m, NbParticles, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R = TransferDataClass.OpenMextract(PathMextract)
    ToolWidget = SpaceView(t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R)
    # ToolWidget = SpaceView(t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R)
    # ToolWidget = DiagramAE(t_m, NbBodies_m, a_m, e_m)
    ToolWidget.show() 
    app.exec() # Application execution
