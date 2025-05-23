#! /Users/lacquema/Astroide.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os
import numpy as np
from numpy import cos, sin, exp, log, log10, linspace, max, loadtxt, transpose, histogram, histogram2d, arcsinh, array, float64, int8
from math import pi, sqrt
from random import random

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QCheckBox

# My packages
from WindowParam import WindowParamClass
from WindowPlot import WindowPlot
from Parameters import *
from TransferData import TransferDataClass
from Resonance import ResClass
from Curve import CurveClass

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
        self.WindowPlot.WidgetParam.BtnRefresh.clicked.connect(self.Refresh_active_plots)

        # Widget container
        self.setLayout(self.Layout) # GeneralToolClass is directly the widget container

    # Open the plot window when the Plot button is clicked
    def Toggle_WindowPlot(self):
        self.WindowPlot.show()
        self.BtnPlot.setEnabled(False)
        self.Refresh_active_plots()

    # Refresh all active plot when the refresh button is clicked
    def Refresh_active_plots(self):
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            if WidgetPlot.isVisible():
                WidgetPlot.refresh_plot()
    
    # Change the index of the snapshot
    def Change_IndexSnap(self, value):
        self.IndexSnap = value
        self.Refresh_active_plots()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        app.closeAllWindows()

    # Reset all widgets of the parameters window
    def ResetParams(self):
        for i in reversed(range(2, self.WindowPlot.WidgetParam.Layout.count())): 
            WidgetToRemove = self.WindowPlot.WidgetParam.Layout.itemAt(i).widget()
            self.WindowPlot.WidgetParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        self.InitParams()
        for WidgetPlot in self.WindowPlot.WidgetPlots:
            WidgetPlot.reset_history()
        self.Refresh_active_plots()

    # Compute of the limits of subplot 2d
    def subplot_lim_2d(self, widget_plot, xlim_init=None, ylim_init=None):
        """
        Determine the axis limits for a subplot 2d based on the widget's history.
        """
        if not widget_plot.history:  # If history is empty
            return xlim_init, ylim_init
        else:
            # Retrieve limits from history
            index = widget_plot.history_index
            xlim = widget_plot.history[index].get('xlim', xlim_init)
            ylim = widget_plot.history[index].get('ylim', ylim_init)
            return xlim, ylim
        
    # Compute of the limits of subplot 2d
    def subplot_lim_3d(self, widget_plot, xlim_init=None, ylim_init=None, zlim_init=None):
        """
        Determine the axis limits for a subplot 3d based on the widget's history.
        """
        if not widget_plot.history:  # If history is empty
            return xlim_init, ylim_init, zlim_init
        else:
            # Retrieve limits from history
            index = widget_plot.history_index
            xlim = widget_plot.history[index].get('xlim', xlim_init)
            ylim = widget_plot.history[index].get('ylim', ylim_init)
            zlim = widget_plot.history[index].get('zlim', zlim_init)
            return xlim, ylim, zlim

    def InitParams(self):
        return

    def Plot(self):
        return
    


class SpaceView(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R):
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

        # Parameters initialisation
        self.InitParams()
        self.InitWidgetPlots()

    def InitWidgetPlots(self):
        self.WidgetPlotXY = self.WindowPlot.add_WidgetPlot(self.PlotXY)
        self.WidgetPlotXZ = self.WindowPlot.add_WidgetPlot(self.PlotXZ)
        self.WidgetPlotXYZ = self.WindowPlot.add_WidgetPlot(self.PlotXYZ)
        self.indexViewChanged(self.indexView)
        
    def InitParams(self):
        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter(Title='View :'))

        # Type of view
        self.ViewWidget = ComboBox('View', 'Dimension', ['face-on', 'edge-on', '3D'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ViewWidget)
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.ViewWidget.ComboParam.currentIndexChanged.connect(self.indexViewChanged)

        # X limits
        self.LimDefault = int(round(max(self.a_m[0])*(1+max(self.e_m[0]))))

        # self.Xmin = -LimDefault
        # self.XminWidget = SpinBox('Xmin', 'X minimum [AU]', self.Xmin)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.XminWidget)

        # self.Xmax = +LimDefault
        # self.XmaxWidget = SpinBox('Xmax', 'X maximum [AU]', self.Xmax)
        # self.XminWidget.Layout.addWidget(self.XmaxWidget) 

        # Y limits
        # self.Ymin = -LimDefault
        # self.YminWidget = SpinBox('Ymin', 'Y minimum [AU]', self.Ymin)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.YminWidget)

        # self.Ymax = +LimDefault
        # self.YmaxWidget = SpinBox('Ymax', 'Y maximum [AU]', self.Ymax)
        # self.YminWidget.Layout.addWidget(self.YmaxWidget)

        # Z limits
        # self.Zmin = -LimDefault
        # self.ZminWidget = SpinBox('Zmin', 'Z minimum [AU]', self.Zmin)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.ZminWidget)
        # self.ZminWidget.setEnabled(False)

        # self.Zmax = +LimDefault
        # self.ZmaxWidget = SpinBox('Zmax', 'Z maximum [AU]', self.Zmax)
        # self.ZminWidget.Layout.addWidget(self.ZmaxWidget)
        # self.ZmaxWidget.setEnabled(False)

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
        
        self.Refresh_active_plots()


        # if self.indexView == 0:
        #     self.XminWidget.setEnabled(True)
        #     self.XmaxWidget.setEnabled(True)
        #     self.YminWidget.setEnabled(True)
        #     self.YmaxWidget.setEnabled(True)
        #     self.ZminWidget.setEnabled(False)
        #     self.ZmaxWidget.setEnabled(False)
        #     self.RepresWidget.setEnabled(True)

        # elif self.indexView == 1:
        #     self.XminWidget.setEnabled(True)
        #     self.XmaxWidget.setEnabled(True)
        #     self.YminWidget.setEnabled(False)
        #     self.YmaxWidget.setEnabled(False)
        #     self.ZminWidget.setEnabled(True)
        #     self.ZmaxWidget.setEnabled(True)
        #     self.RepresWidget.setEnabled(True)

        # elif self.indexView == 2:
        #     self.XminWidget.setEnabled(True)
        #     self.XmaxWidget.setEnabled(True)
        #     self.YminWidget.setEnabled(True)
        #     self.YmaxWidget.setEnabled(True)
        #     self.ZminWidget.setEnabled(True)
        #     self.ZmaxWidget.setEnabled(True)
        #     self.RepresWidget.setEnabled(False)

    def indexRepresChanged(self, value):
        self.indexRepres = value

        if self.indexRepres == 0:
            self.NbBinsXWidget.setEnabled(False)
            self.NbBinsYWidget.setEnabled(False)
            self.SizePartWidget.setEnabled(True)

        elif self.indexRepres == 1:
            self.NbBinsXWidget.setEnabled(True)
            self.NbBinsYWidget.setEnabled(True)
            self.SizePartWidget.setEnabled(False)

    def UpdateParams(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        # self.Xmin = self.XminWidget.SpinParam.value()
        # self.Xmax = self.XmaxWidget.SpinParam.value()
        # self.Ymin = self.YminWidget.SpinParam.value()
        # self.Ymax = self.YmaxWidget.SpinParam.value()
        # self.Zmin = self.ZminWidget.SpinParam.value()
        # self.Zmax = self.ZmaxWidget.SpinParam.value()
        # self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        self.SizePart = self.SizePartWidget.SpinParam.value()
        self.NbBinsX = self.NbBinsXWidget.SpinParam.value()
        self.NbBinsY = self.NbBinsYWidget.SpinParam.value()

    def general_plot(self):

        # Bodies positions
        self.Xbf, self.Ybf, self.Zbf = [], [], []
        self.Xbm, self.Ybm, self.Zbm = [], [], []
        for k in range(self.NbBodies_m[self.IndexSnap]):
            # xbf, ybf, zbf = self.Ellipse(af[k][isave], ef[k][isave], i[k][isave]*pi/180, W[k][isave]*pi/180, w[k][isave]*pi/180)
            xbm, ybm, zbm = self.Ellipse2(self.a_m[self.IndexSnap][k], self.e_m[self.IndexSnap][k], self.Ex[self.IndexSnap][k], self.Ey[self.IndexSnap][k], self.Ez[self.IndexSnap][k], self.Epx[self.IndexSnap][k], self.Epy[self.IndexSnap][k], self.Epz[self.IndexSnap][k])
            # self.Xbf.append(xbf)
            # self.Ybf.append(ybf)
            # self.Zbf.append(zbf)
            self.Xbm.append(xbm)
            self.Ybm.append(ybm)
            self.Zbm.append(zbm)

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.WidgetPlots[0].Canvas.draw()
            return
        
        # Specific variables
        self.NbBodies = self.NbBodies_m[self.IndexSnap]
        self.t = self.t_m[self.IndexSnap]/10**6

    def PlotXY(self):

        # General plot
        self.general_plot()

        # add subplot
        self.SubplotXY = self.WidgetPlotXY.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')

        # X, Y limits
        (Xmin, Xmax), (Ymin, Ymax) = self.subplot_lim_2d(self.WidgetPlotXY, xlim_init=[-self.LimDefault, self.LimDefault], ylim_init=[-self.LimDefault, self.LimDefault])

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
            hist, xedges, yedges = histogram2d(
                self.X[self.IndexSnap][self.NbBodies-1:], 
                self.Y[self.IndexSnap][self.NbBodies-1:], 
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
        self.SubplotXY.set_xlim(Xmin, Xmax)
        self.SubplotXY.set_ylabel('Y [AU]')
        self.SubplotXY.set_ylim(Ymin, Ymax)


    def PlotXZ(self):

        # General plot
        self.general_plot()

        # add subplot
        self.SubplotXZ = self.WidgetPlotXZ.Canvas.fig.add_subplot(111, aspect='equal', label='Main plot')

        # X, Z limits
        (Xmin, Xmax), (Zmin, Zmax) = self.subplot_lim_2d(self.WidgetPlotXZ, xlim_init=[-self.LimDefault, self.LimDefault], ylim_init=[-self.LimDefault, self.LimDefault])

        # Plot
        for k in range(self.NbBodies):
            if self.CheckOrbits.CheckParam.isChecked():
                self.SubplotXZ.plot(self.Xbm[k], self.Zbm[k], color=self.colorList[k], linestyle='--', label="Orbit of "+str(k+1))
            if self.CheckBodies.CheckParam.isChecked():
                self.SubplotXZ.plot(self.X[self.IndexSnap][k], self.Z[self.IndexSnap][k], marker='.', markersize=self.SizeBodies, color=self.colorList[k], label="Marker of "+str(k+1))

        if self.indexRepres == 0:
            self.SubplotXZ.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Z[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0, label='Particles')
        
        elif self.indexRepres == 1:
            hist, xedges, yedges = histogram2d(
                self.X[self.IndexSnap][self.NbBodies-1:], 
                self.Z[self.IndexSnap][self.NbBodies-1:], 
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
        self.SubplotXZ.set_xlim(Xmin, Xmax)
        self.SubplotXZ.set_ylabel('Z [AU]')
        self.SubplotXZ.set_ylim(Zmin, Zmax)
        

    def PlotXYZ(self):

        # General plot
        self.general_plot()

        # Add subplot
        self.SubplotXYZ = self.WidgetPlotXYZ.Canvas.fig.add_subplot(111, aspect='equal', projection='3d', label='Main plot')

        # X, Y, Z limits
        (Xmin, Xmax), (Ymin, Ymax), (Zmin, Zmax) = self.subplot_lim_3d(self.WidgetPlotXYZ, xlim_init=[-self.LimDefault, self.LimDefault], ylim_init=[-self.LimDefault, self.LimDefault], zlim_init=[-self.LimDefault, self.LimDefault])

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

        # Ensure axis limits are finite
        self.SubplotXYZ.set_xlim(Xmin, Xmax)
        self.SubplotXYZ.set_ylim(Ymin, Ymax)
        self.SubplotXYZ.set_zlim(Zmin, Zmax)

        # # Set aspect ratio to 'auto' to avoid divide-by-zero errors
        # self.SubplotXYZ.set_box_aspect([1, 1, 1])  # Equal aspect ratio for all axes

        # Plot features
        self.SubplotXYZ.set_xlabel('X [AU]')
        self.SubplotXYZ.set_ylabel('Y [AU]')
        self.SubplotXYZ.set_zlabel('Z [AU]')

        # Time
        self.SubplotXYZ.text2D(1.01, 1.01, 't='+str(round(self.t, 1))+' Myr', 
                               horizontalalignment='right', verticalalignment='top', 
                               transform=self.SubplotXYZ.transAxes, label='Time')

        # Add legend
    def Ellipse2(self, a, e, Ex, Ey, Ez, Epx, Epy, Epz):

        E = linspace(-pi, pi, 100)
        x = Ex*a*(cos(E)-e) + Epx*a*sqrt(1-e**2)*sin(E)
        y = Ey*a*(cos(E)-e) + Epy*a*sqrt(1-e**2)*sin(E)
        z = Ez*a*(cos(E)-e) + Epz*a*sqrt(1-e**2)*sin(E)

        return x, y, z


class RadProfile(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, X, Y, Z, R):
        super().__init__('Radial profile', "Particules' integrated radial profile")

        # Data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.R = R
        self.X = X
        self.Y = Y
        self.Z = Z 

        # Beta Pic other curves
        self.profileAug = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Augereau_profile.dat', dtype = float))
        self.profileNE = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Dent_profile_NE.dat', dtype = float))
        self.profileSW = transpose(loadtxt(self.DirPath+'/OtherCurves/bpic_Dent_profile_SW.dat', dtype = float))


        # Parameters initialisation
        self.InitParams()

        # Plot initialisation
        events_to_reset_history = [self.Norm.ComboParam.currentIndexChanged, self.Ordinate.ComboParam.currentIndexChanged]
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, events_to_reset_history)

    def InitParams(self):

        # # R limits
        self.Rmin = 0
        # self.RminWidget = SpinBox('Rmin', 'Radii minimum [AU]', self.Rmin)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.RminWidget)

        # # print(self.R)
        self.Rmax = round(max(self.R[0]))
        # self.RmaxWidget = SpinBox('Rmax', 'Radii maximum [AU]', self.Rmax)
        # self.RminWidget.Layout.addWidget(self.RmaxWidget) 

        # self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Histogram
        self.Ordinate = ComboBox('Ordinate', 'Choice of ordinate', ['Number of particules', 'Surface density'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.Ordinate)

        self.Norm = ComboBox('Normalisation', 'Choice of the normalisation', ['None', 'MaxEqOne', 'SumEqOne'])
        self.Norm.ComboParam.setCurrentIndex(0)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.Norm)

        self.NbBins = 100
        self.NbBinsWidget = SpinBox('Bining', 'Number of bins', self.NbBins, 1)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.NbBinsWidget)

        # Checkbox pour l'interpolation
        self.CheckInterpolation = CheckBox('Interpolation')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckInterpolation)


        # self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

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

        self.CheckAugWidget = QCheckBox('Aug 2001')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckAugWidget)
        self.CheckAugWidget.setEnabled(self.CheckCurves.isChecked())

        self.CheckDentWidget = QCheckBox('Dent 2014')
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckDentWidget)
        self.CheckDentWidget.setEnabled(self.CheckCurves.isChecked())

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
        self.CheckAugWidget.setEnabled(CurveState)
        self.CheckDentWidget.setEnabled(CurveState)
        self.ButAddCurve.setEnabled(CurveState)
        for i in range(len(self.CurveWidgets)):
            self.CurveWidgets[i].setEnabled(CurveState)
    
    def UpdateParams(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        # self.Rmin = self.RminWidget.SpinParam.value()
        # self.Rmax = self.RmaxWidget.SpinParam.value()
        self.NbBins = self.NbBinsWidget.SpinParam.value()
        # self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        if self.CheckCurves.isChecked():
            self.CurvePaths = []
            self.CurveLabels = []
            for i in range(len(self.CurveWidgets)):
                self.CurvePaths.append(self.CurveWidgets[i].PathWidget.EditParam.text())
                self.CurveLabels.append(self.CurveWidgets[i].LabelWidget.EditParam.text())
                if self.CurveLabels[i] == '': self.CurveLabels[i] = self.CurvePaths[i].split('/')[-1]

    def Rebin(self, data, resolution):    
        """
        Rebin the data to a specified resolution.
        """
        smooth_data, bins_edges = histogram(data, bins=(np.max(data)-np.min(data))/resolution)
        return smooth_data 

    # Plot
    def Plot(self):

        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # X, Z limits
        (Rmin, Rmax), _ = self.subplot_lim_2d(self.WidgetPlot, xlim_init=[-self.Rmin, self.Rmax])

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WidgetPlot.Canvas.draw()
            return
        
        # Specific variables
        self.t = self.t_m[self.IndexSnap]/10**6

        # Plot with current parameters
        if self.CheckBodies.CheckParam.isChecked():
            for k in range(self.NbBodies_m[self.IndexSnap]):
                self.Subplot.plot(self.a_m[self.IndexSnap][k], 0, color=self.colorList[k], marker='.', markersize=self.SizeBodies)

        # Histogram
        # X_smoothed = gaussian_filter(self.X[self.IndexSnap], sigma=1)
        # Y_smoothed = gaussian_filter(self.Y[self.IndexSnap], sigma=1)
        # Z_smoothed = gaussian_filter(self.Z[self.IndexSnap], sigma=1)
        # self.R = np.sqrt(self.X[self.IndexSnap]**2+self.Y[self.IndexSnap]**2+self.Z[self.IndexSnap]**2)
        histCount, histX = histogram([x for x in self.R[self.IndexSnap]  if Rmin<x<Rmax], bins=self.NbBins)

        # Surface density computatiom
        if self.Ordinate.ComboParam.currentIndex()==1: histCount = histCount/(2*pi*histX[:-1]) # surface density

        # Normalisation computation
        self.NormDiv = 1
        if self.Norm.ComboParam.currentIndex()==1: self.NormDiv = max(histCount) # normalisation of max equal one
        elif self.Norm.ComboParam.currentIndex()==2: self.NormDiv = sum(histCount) # normalisation of sum equal one
        histCount = histCount/self.NormDiv

        # Interpolation si la checkbox est cochée
        if self.CheckInterpolation.CheckParam.isChecked():
            try:
                # Utilisation des milieux des bins pour l'interpolation
                histX_mid = (histX[:-1] + histX[1:]) / 2  # Calcul des milieux des bins
                interp_func = interp1d(histX_mid, histCount, kind='cubic', fill_value="extrapolate")
                
                # Génération de points interpolés
                histX_interp = linspace(histX_mid[0], histX_mid[-1], 1000)  # Plus de points pour une courbe lisse
                histCount_interp = interp_func(histX_interp)
                
                # Tracé de la courbe interpolée
                self.Subplot.plot(histX_interp, histCount_interp, label='Interpolation', linewidth=1, color='green')
            except Exception as e:
                print(f"Interpolation failed: {e}")

        # Stairs
        self.Subplot.stairs(histCount, histX, label='Simulation', linewidth=1, color='black')
        
        # Other curves 
        if self.CheckAugWidget.isChecked(): self.Subplot.plot(self.profileAug[0], self.profileAug[1], color='blue', linestyle='dashed', linewidth=0.5, label='Aug+2001')
        if self.CheckDentWidget.isChecked(): self.Subplot.plot((self.profileNE[0]+self.profileSW[0])/2, (self.profileNE[1]+self.profileSW[1])/2, color='red', linestyle='dashed', linewidth=0.5, label='Dent+2014')
        
        if self.CheckCurves.isChecked(): 
            try:
                for i in range(len(self.CurveWidgets)):
                    self.Curve = transpose(np.loadtxt(self.CurvePaths[i], dtype=float, delimiter=',' if ',' in open(self.CurvePaths[i]).readline() else '\t'))
                    self.Subplot.plot(self.Curve[0], self.Curve[1], linestyle='dashed', linewidth=1, label=self.CurveLabels[i])
            except:
                print('Wrong files')
                print('Path is necessary to files in format:')
                print('Column 1 : abscissa')
                print('Column 2 : ordinate')

        # Time
        self.Subplot.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.Subplot.transAxes)

        # Plot features
        if self.CheckCurves.isChecked(): self.Subplot.legend()
        self.Subplot.set_xlabel('Radius [AU]')
        self.Subplot.set_xlim(Rmin, Rmax)
        if self.Ordinate.ComboParam.currentIndex()==0: self.Subplot.set_ylabel('Number of particules')
        elif self.Ordinate.ComboParam.currentIndex()==1: self.Subplot.set_ylabel('Surface density [arbitrary unit]')
        self.Subplot.set_ylim(0, None)
        # self.Subplot.set_ylim(0, 1)
        
        # Update canvas
        # self.WindowPlot.WidgetPlots.Canvas.draw()



class DiagramAE(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m):
        super().__init__('Diagram a=f(e)', "Bodies' evolution of semi-major axis as a function of excentricity")

        # General data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.e_m = e_m

        # Parameters initialisation
        self.InitParams()

        # Plot initialisation
        events_to_reset_history = []
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, events_to_reset_history)

    def InitParams(self):

        # Limits
        self.Rmin = 0
        # self.RminWidget = SpinBox('Rmin', 'Minimum radius [AU]', self.Rmin, 0, None)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.RminWidget)

        self.Rmax = round(max(self.a_m[0]))
        # self.RmaxWidget = SpinBox('Rmax', 'Maximum radius [AU]', self.Rmax, 0, None)
        # self.RminWidget.Layout.addWidget(self.RmaxWidget) 

        self.Emin = 0
        # self.EminWidget = DoubleSpinBox('Emin', 'Minimum eccentricity', self.Emin, 0, None)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.EminWidget)

        self.Emax = 1
        # self.EmaxWidget = DoubleSpinBox('Emax', 'Maximum eccentricity', self.Emax, 0, None)
        # self.WindowPlot.WidgetParam.Layout.addWidget(self.EmaxWidget) 

        # self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

        # Bodies' positions
        self.CheckBodies = CheckBox("Bodies' position showing")
        self.WindowPlot.WidgetParam.Layout.addWidget(self.CheckBodies, alignment=Qt.AlignmentFlag.AlignCenter)

        self.SizeBodies = 15
        self.SizeBodiesWidget = SpinBox('Size of bodies', 'Plot size of bodies', self.SizeBodies)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.SizeBodiesWidget)
        
        self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked())
        self.CheckBodies.CheckParam.stateChanged.connect(lambda: self.SizeBodiesWidget.setEnabled(self.CheckBodies.CheckParam.isChecked()))

        self.WindowPlot.WidgetParam.Layout.addWidget(Delimiter())

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
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        # self.Rmin = self.RminWidget.SpinParam.value()
        # self.Rmax = self.RmaxWidget.SpinParam.value()
        # self.Emin = self.EminWidget.SpinParam.value()
        # self.Emax = self.EmaxWidget.SpinParam.value()
        self.SizePart = self.SizePartWidget.SpinParam.value()
        self.SizeBodies = self.SizeBodiesWidget.SpinParam.value()
        if self.CheckRes.isChecked():
            self.nRefValues = []
            self.PResValues = []
            self.PRefValues = []
            for i in range(len(self.ResWidgets)):
                self.nRefValues.append(int(self.ResWidgets[i].nRefWidget.SpinParam.value()))
                self.PResValues.append(int(self.ResWidgets[i].PResWidget.SpinParam.value()))
                self.PRefValues.append(int(self.ResWidgets[i].PRefWidget.SpinParam.value()))
        
    # Plot
    def Plot(self):
        
        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WidgetPlot.Canvas.draw()
            return
    
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
        
        # Plot with current parameters
        if self.CheckBodies.CheckParam.isChecked():
            for k in range(self.NbBodies_m[self.IndexSnap]):
                self.Subplot.plot(self.a[k], self.e[k], color=self.colorList[k], marker='.', markersize=self.SizeBodies)
        self.Subplot.scatter(self.a[self.NbBodies_m[self.IndexSnap]:], self.e[self.NbBodies_m[self.IndexSnap]:], s=self.SizePart, c='black', linewidths=0)

        if self.CheckRes.isChecked():
            for i in range(len(self.ResWidgets)):
                colorRef = self.colorList[self.nRefValues[i]]
                self.Subplot.axvline(self.aResValues[i], linewidth=1, color=colorRef, linestyle='--')
                # self.Subplot.annotate(str(self.PResValues[i])+':'+str(self.PRefValues[i]), xy=(self.aResValues[i], 0.9), bbox=dict(boxstyle='round', facecolor='white', edgecolor=colorRef))
                self.Subplot.text(self.aResValues[i], 0.9, ' '+str(self.PResValues[i])+':'+str(self.PRefValues[i])+' ', rotation=90, color=colorRef, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=colorRef), fontsize=7, va='center', ha='center')

        # Time
        self.Subplot.text(x=0.99, y=1.01, s='t='+str(round(self.t, 1))+' Myr', horizontalalignment='right', verticalalignment='bottom', transform=self.Subplot.transAxes)

        # Plot features
        # self.Subplot.set_title('t='+str(round(self.t, 1))+' Myr')
        self.Subplot.set_xlabel('Semi-major axis [AU]')
        self.Subplot.set_xlim(self.Rmin, self.Rmax)
        self.Subplot.set_ylabel('Eccentricity')
        self.Subplot.set_ylim(self.Emin, self.Emax)
            
        # Update canvas
        # self.WindowPlot.WidgetPlots.Canvas.draw()



class DiagramTY(GeneralToolClass):
    def __init__(self, NbBodies_f, t_f, a_f, e_f, i, W, w, M):
        super().__init__('Diagram y=f(t)', "Planets' evolution of an orbital parameter as a function of time")

        # General data
        self.NbBodies_f = NbBodies_f
        self.t_f = t_f
        self.a_f = a_f
        self.e_f = e_f
        self.i = i
        self.W = W
        self.w = w
        self.M = M

        # Parameters initialisation
        self.InitParams()

        # Plot initialisation
        events_to_reset_history = [self.ParamOrbitWidget.ComboParam.currentIndexChanged]
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, events_to_reset_history)

    def InitParams(self):
        # self.CheckTitle.setEnabled(False)

        # Orbit number
        self.nOrbit = 1
        self.nOrbitWidget = SpinBox("Bodie's number",'Number of the body counting outwards', ParamMin=1, ParamMax=self.NbBodies_f-1)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.nOrbitWidget)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['a','e','i','W','w','M'])
        self.WindowPlot.WidgetParam.Layout.addWidget(self.ParamOrbitWidget)

        # Time limits
        self.Tmin = 0
        self.TminWidget = SpinBox('Tmin', 'Time minimum [yr]', self.Tmin)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.TminWidget)

        self.Tmax = round(max(self.t_f))
        self.TmaxWidget = SpinBox('Tmax', 'Time maximum [yr]', self.Tmax)
        self.TminWidget.Layout.addWidget(self.TmaxWidget) 

    def UpdateParams(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        self.nOrbit = self.nOrbitWidget.SpinParam.value()
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.Tmin = self.TminWidget.SpinParam.value()
        self.Tmax = self.TmaxWidget.SpinParam.value()

    # Plot
    def Plot(self):
        
        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WidgetPlot.Canvas.draw()
            return
        
        # Specific data
        t = [self.t_f[self.nOrbit],'Time [yr]']
        a = [self.a_f[self.nOrbit],'Semi-major axis [AU]']
        e = [self.e_f[self.nOrbit],'Eccentricity']
        i = [self.i[self.nOrbit],'Inclinaison [deg]']
        W = [self.W[self.nOrbit],'Longitude of ascending node [deg]']
        w = [self.w[self.nOrbit],'Argument of periastron [deg]']
        M = [self.M[self.nOrbit],'Initial mean longitude']
        
        # Plot with current parameters
        self.EvalParamOrbit = eval(self.ParamOrbit)
        self.Subplot.plot(t[0], self.EvalParamOrbit[0], linewidth=0.5)

        # Plot features
        self.Subplot.set_xlabel('Time [yr]')
        self.Subplot.set_xlim(self.Tmin, self.Tmax)
        self.Subplot.set_ylabel(self.EvalParamOrbit[1])
        
        # # Update canvas
        # self.WindowPlot.WidgetPlots.Canvas.draw()



class DiagramXY(GeneralToolClass):
    def __init__(self, NbBodies_f, t_f, a_f, e_f, i, W, w, M):
        super().__init__('Diagram y=f(x)', "Planets' evolution of an orbital parameter Y as a function of an other X")

        # General data
        self.NbBodies_f = NbBodies_f
        self.t_f = t_f
        self.a_f = a_f
        self.e_f = e_f
        self.i = i
        self.W = W
        self.w = w
        self.M = M

        # Parameters initialisation
        self.InitParams()

        # Plot initialisation
        events_to_reset_history = []
        self.WidgetPlot = self.WindowPlot.add_WidgetPlot(self.Plot, events_to_reset_history)

    def InitParams(self):
        # self.CheckTitle.setEnabled(False)

        # Ordinate quantity formula
        self.YFormula = ''
        self.YFormulaWidget = LineEdit('Y formula','Ordinate quantity by combining t, a, e, i, w, W, [n] where n is the number of the body counting outwards, and math fonctions', self.YFormula)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YFormulaWidget)
        self.YFormulaWidget.setMinimumWidth(300)

        # Ordinate quantity label
        self.YLabel = ''
        self.YLabelWidget = LineEdit('Y label', 'Ordinate label', self.YLabel)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.YLabelWidget)
        # self.YLabelWidget.setEnabled(self.CheckYLabel.CheckParam.isChecked())

        # Abscissa quantity formula
        self.XFormula = '' 
        self.XFormulaWidget = LineEdit('X formula','Abscissa quantity by combining t, a, e, i, w, W, [n] where n is the number of the body counting outwards, and math fonctions', self.XFormula)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XFormulaWidget)
        # self.XFormulaWidget.EditParam.textEdited.connect(self.Validator)
        self.XFormulaWidget.setMinimumWidth(300)

        # Abscissa quantity label
        self.XLabel = ''
        self.XLabelWidget = LineEdit('X label', 'Abscissa label', self.XLabel)
        self.WindowPlot.WidgetParam.Layout.addWidget(self.XLabelWidget)
        # self.XLabelWidget.setEnabled(self.CheckXLabel.CheckParam.isChecked())

        # Connections
        # self.CheckXLabel.CheckParam.stateChanged.connect(lambda: self.XLabelWidget.setEnabled(self.CheckXLabel.CheckParam.isChecked()))
        # self.CheckYLabel.CheckParam.stateChanged.connect(lambda: self.YLabelWidget.setEnabled(self.CheckYLabel.CheckParam.isChecked()))

    # def Validator(self, text):
    #     for pos in range(1, len(text)): 
    #         if text[pos-1]=='[' and not 1<=text[pos]<=self.NbBodies_f-1:
    #             self.XFormulaWidget.EditParam.setText(self.TextOld)


    def UpdateParameters(self):
        # self.SizeLabels = self.SizeLabelsWidget.SpinParam.value()
        self.YFormula = self.YFormulaWidget.EditParam.text()
        self.XFormula = self.XFormulaWidget.EditParam.text()
        self.XLabel = self.XLabelWidget.EditParam.text()

        # # Mofification of orbit count on X-formula
        # orbit_num = []
        # orbit_num_pos = []
        # for pos in range(1, len(self.XFormula)): 
        #     if self.XFormula[pos-1]=='[':
        #         orbit_num.append(int(self.XFormula[pos]))
        #         orbit_num_pos.append(pos)
        # orbit_num = self.NbBodies_f-1-array(orbit_num)
        # for j in range(len(orbit_num)): 
        #     self.XFormula = list(self.XFormula)
        #     self.XFormula[orbit_num_pos[j]] = str(orbit_num[j])
        #     self.XFormula = ''.join(self.XFormula)

        # # Mofification of orbit count on Y-formula
        # orbit_num = []
        # orbit_num_pos = []
        # for pos in range(1, len(self.YFormula)): 
        #     if self.YFormula[pos-1]=='[':
        #         orbit_num.append(int(self.YFormula[pos]))
        #         orbit_num_pos.append(pos)
        # orbit_num = self.NbBodies_f-1-array(orbit_num)
        # for j in range(len(orbit_num)): 
        #     self.YFormula = list(self.YFormula)
        #     self.YFormula[orbit_num_pos[j]] = str(orbit_num[j])
        #     self.YFormula = ''.join(self.YFormula)

        # Modification of labels
        if len(self.YLabelWidget.EditParam.text()) != 0: self.YLabel = self.YLabelWidget.EditParam.text()
        else: self.YLabel = self.YFormula
        if len(self.XLabelWidget.EditParam.text()) != 0: self.XLabel = self.XLabelWidget.EditParam.text()
        else: self.XLabel = self.XFormula
                

    # Plot
    def Plot(self):

        # Plot initialisation
        self.Subplot = self.WidgetPlot.Canvas.fig.add_subplot(111)

        # Update of parameters
        try:
            self.UpdateParameters()
        except:
            print('Wrong Parameters')
            self.WidgetPlot.Canvas.draw()
            return

        # Specific variables
        t = self.t_f
        a = self.a_f
        e = self.e_f
        i = self.i
        W = self.W
        w = self.w
        M = self.M

        # Plot with current parameters
        if self.XFormula != '' and self.YFormula != '':
            try:
                self.Y = eval(self.YFormula)
                self.X = eval(self.XFormula)
                self.Subplot.plot(self.X, self.Y, linewidth=0.5)
            except:
                print('Wrong formula ! Use t, a, e, i, w, W, [n] where n is the number of orbit, and math fonctions')
        else:
            print('Wrong formula ! Use t, a, e, i, w, W, [n] where n is the number of orbit, and math fonctions')

        # Plot features
        self.Subplot.set_xlabel(self.XLabel)
        self.Subplot.set_ylabel(self.YLabel)
        
        # # Update canvas
        # self.WidgetPlot.Canvas.draw()






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
