#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
from numpy import cos, sin, exp, log, log10, linspace, max, loadtxt, transpose, histogram, histogram2d, arcsinh
from math import pi, sqrt
from random import random

# PyQt packages
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QCheckBox

# My packages
from WindowParam import WindowParamClass
from WindowPlot import WindowPlotClass
from Parameters import *
from TransferData import TransferDataClass
from Resonance import ResClass
from Curve import CurveClass

from mpl_toolkits.axes_grid1 import make_axes_locatable


### --- Tools Generating --- ###

class GeneralToolClass(QWidget):

    def __init__(self, ToolName, ToolStatus):
        super().__init__()

        self.IndexSnap = 0
        self.colorList = ['black', 'blue', 'red', 'green', 'orange', 'pink']

        # Layout
        Layout = QHBoxLayout()

        # Label
        LblTool = QLabel("{} :".format(ToolName))
        LblTool.setStatusTip(ToolStatus)
        Layout.addWidget(LblTool)

        # Parameters button
        self.BtnParam = QPushButton('Parameters')
        self.BtnParam.clicked.connect(self.Toggle_WindowParam)
        Layout.addWidget(self.BtnParam)

        # Initialisation of parameters windows
        self.WindowParam = WindowParamClass(ToolName)
        self.WindowParam.BtnReset.clicked.connect(self.ResetParams)
        self.WindowParam.SignalCloseWindowParam.connect(lambda: self.BtnParam.setEnabled(True)) # reception of the closeEvent of the parameter window and set enabled the associed button
        self.WindowParam.resize(self.WindowParam.minimumSize())

        # Plot button
        self.BtnPlot = QPushButton('Plot')
        self.BtnPlot.clicked.connect(self.Toggle_WindowPlot)
        Layout.addWidget(self.BtnPlot)

        # Initialisation of plot windows
        self.WindowPlot = WindowPlotClass(ToolName)
        self.WindowPlot.SignalCloseWindowPlot.connect(lambda: self.BtnPlot.setEnabled(True)) # reception of the closeEvent of the plot window and set enabled the associed button

        # Widget container
        self.setLayout(Layout) # GeneralToolClass is directly the widget container


    # Open the parameters window when the parameters button is clicked
    def Toggle_WindowParam(self):
            # self.WindowParam.move(self.x()-500, self.pos().y())
            self.WindowParam.show()
            self.BtnParam.setEnabled(False)

    # Open the plot window when the Plot button is clicked
    def Toggle_WindowPlot(self):
            self.Plot()
            self.WindowPlot.show()
            self.BtnPlot.setEnabled(False)

    # Refresh all active plot when the refresh button is clicked
    def Refresh_ActivePlots(self):
        if self.WindowPlot.isVisible():
            self.Plot()
            
    def Change_IndexSnap(self, value):
        self.IndexSnap = value
        self.Refresh_ActivePlots()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        app.closeAllWindows()

    # Reset all widgets of the parameters window
    def ResetParams(self):
        for i in reversed(range(1, self.WindowParam.Layout.count())): 
            WidgetToRemove = self.WindowParam.Layout.itemAt(i).widget()
            self.WindowParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        self.InitParams()

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
        
        # # Plots initialisation
        # self.Subplot3D = self.WindowPlot.Canvas.fig.add_subplot(111, projection='3d', aspect='equal')
        # self.Subplot2D = self.WindowPlot.Canvas.fig.add_subplot(111, aspect='equal')

        # Parameters initialisation
        self.InitParams()


    def InitParams(self):
        # Type of view
        self.ViewWidget = ComboBox('View', 'Dimension', ['face-on', 'edge-on', '3D'])
        self.WindowParam.Layout.addWidget(self.ViewWidget)
        self.indexView = self.ViewWidget.ComboParam.currentIndex()
        self.ViewWidget.ComboParam. currentIndexChanged.connect(self.indexViewChanged)

        # Planets' positions
        self.BodPosWidget = QCheckBox("Planets' position showing")
        self.WindowParam.Layout.addWidget(self.BodPosWidget)

        # X limits
        LimDefault = int(round(max(self.a_m[0])*(1+max(self.e_m[0]))))

        self.Xmin = -LimDefault
        self.XminWidget = LineEditInt('Xmin', 'X minimum(AU)', self.Xmin)
        self.WindowParam.Layout.addWidget(self.XminWidget)

        self.Xmax = +LimDefault
        self.XmaxWidget = LineEditInt('Xmax', 'X maximum (AU)', self.Xmax)
        self.XminWidget.Layout.addWidget(self.XmaxWidget) 

        # Y limits
        self.Ymin = -LimDefault
        self.YminWidget = LineEditInt('Ymin', 'Y minimum (AU)', self.Ymin)
        self.WindowParam.Layout.addWidget(self.YminWidget)

        self.Ymax = +LimDefault
        self.YmaxWidget = LineEditInt('Ymax', 'Y maximum (AU)', self.Ymax)
        self.YminWidget.Layout.addWidget(self.YmaxWidget)

        # Z limits
        self.Zmin = -LimDefault
        self.ZminWidget = LineEditInt('Zmin', 'Z minimum (AU)', self.Zmin)
        self.WindowParam.Layout.addWidget(self.ZminWidget)
        self.ZminWidget.setEnabled(False)

        self.Zmax = +LimDefault
        self.ZmaxWidget = LineEditInt('Zmax', 'Z maximum (AU)', self.Zmax)
        self.ZminWidget.Layout.addWidget(self.ZmaxWidget)
        self.ZmaxWidget.setEnabled(False)

        # Type of representation
        self.RepresWidget = ComboBox('Representation', 'Type of representation', ['scatter', 'density'])
        self.WindowParam.Layout.addWidget(self.RepresWidget)
        self.indexRepres = self.RepresWidget.ComboParam.currentIndex()
        self.RepresWidget.ComboParam.currentIndexChanged.connect(self.indexRepresChanged)

        self.SizePart = 0.3
        self.SizePartWidget = LineEditFloatPos('Size of particules', 'Plot size of particles', self.SizePart)
        self.WindowParam.Layout.addWidget(self.SizePartWidget)

        self.NbBinsX = 300
        self.NbBinsXWidget = LineEditIntPos('X Bining', 'Number of bins in the x-axis', self.NbBinsX)
        self.WindowParam.Layout.addWidget(self.NbBinsXWidget)
        self.NbBinsXWidget.setEnabled(False)

        self.NbBinsY = 300
        self.NbBinsYWidget = LineEditIntPos('Y Bining', 'Number of bins in the y-axis', self.NbBinsY)
        self.NbBinsXWidget.Layout.addWidget(self.NbBinsYWidget)
        self.NbBinsYWidget.setEnabled(False)


    def indexViewChanged(self, value):
        self.indexView = value

        if self.indexView == 0:
            self.XminWidget.setEnabled(True)
            self.XmaxWidget.setEnabled(True)
            self.YminWidget.setEnabled(True)
            self.YmaxWidget.setEnabled(True)
            self.ZminWidget.setEnabled(False)
            self.ZmaxWidget.setEnabled(False)
            self.RepresWidget.setEnabled(True)

        elif self.indexView == 1:
            self.XminWidget.setEnabled(True)
            self.XmaxWidget.setEnabled(True)
            self.YminWidget.setEnabled(False)
            self.YmaxWidget.setEnabled(False)
            self.ZminWidget.setEnabled(True)
            self.ZmaxWidget.setEnabled(True)
            self.RepresWidget.setEnabled(True)

        elif self.indexView == 2:
            self.XminWidget.setEnabled(True)
            self.XmaxWidget.setEnabled(True)
            self.YminWidget.setEnabled(True)
            self.YmaxWidget.setEnabled(True)
            self.ZminWidget.setEnabled(True)
            self.ZmaxWidget.setEnabled(True)
            self.RepresWidget.setEnabled(False)

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
        self.Xmin = int(self.XminWidget.EditParam.text())
        self.Xmax = int(self.XmaxWidget.EditParam.text())
        self.Ymin = int(self.YminWidget.EditParam.text())
        self.Ymax = int(self.YmaxWidget.EditParam.text())
        self.Zmin = int(self.ZminWidget.EditParam.text())
        self.Zmax = int(self.ZmaxWidget.EditParam.text())
        self.SizePart = float(self.SizePartWidget.EditParam.text().replace(',', '.'))
        self.NbBinsX = int(self.NbBinsXWidget.EditParam.text())
        self.NbBinsY = int(self.NbBinsYWidget.EditParam.text())


    def Plot(self):

        # Clear axis
        for i in range(len(self.WindowPlot.Canvas.fig.axes)):
            self.WindowPlot.Canvas.fig.delaxes(self.WindowPlot.Canvas.fig.axes[0])

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
            self.WindowPlot.Canvas.draw()
            return
        
        # Specific variables
        self.NbBodies = self.NbBodies_m[self.IndexSnap]
        self.t = self.t_m[self.IndexSnap]/10**6


        if self.indexView == 0:

            # Initialisation of axis
            self.Subplot2D = self.WindowPlot.Canvas.fig.add_subplot(111, aspect='equal')

            # Plot
            for k in range(self.NbBodies):
                # self.Subplot2D.plot(self.Xbf[k], self.Ybf[k])
                self.Subplot2D.plot(self.Xbm[k], self.Ybm[k], color=self.colorList[k], linestyle='--')
                if self.BodPosWidget.isChecked():
                    self.Subplot2D.plot(self.X[self.IndexSnap][k], self.Y[self.IndexSnap][k], marker='.', markersize=10, color=self.colorList[k])

            if self.indexRepres == 0:
                self.Subplot2D.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Y[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0)
            elif self.indexRepres == 1:
                hist = self.Subplot2D.hist2d(self.X[self.IndexSnap][self.NbBodies-1:], self.Y[self.IndexSnap][self.NbBodies-1:], bins=[self.NbBinsX, self.NbBinsY], range=[[self.Xmin, self.Xmax],[self.Ymin, self.Ymax]], cmap='magma')
                ColorbarAx = make_axes_locatable(self.Subplot2D).append_axes('right', size='5%', pad=0.1)
                self.WindowPlot.Canvas.fig.colorbar(hist[3], ColorbarAx)

            # Plot features
            self.Subplot2D.set_title('t='+str(round(self.t, 1))+' Myr')
            self.Subplot2D.set_xlabel('X (AU)')
            self.Subplot2D.set_xlim(self.Xmin, self.Xmax)
            self.Subplot2D.set_ylabel('Y (AU)')
            self.Subplot2D.set_ylim(self.Ymin, self.Ymax)



        elif self.indexView == 1:

            # Initialisation of axis
            self.Subplot2D = self.WindowPlot.Canvas.fig.add_subplot(111, aspect='equal')

            # Plot
            for k in range(self.NbBodies):
                self.Subplot2D.plot(self.Xbm[k], self.Zbm[k], color=self.colorList[k], linestyle='--')
                if self.BodPosWidget.isChecked():
                    self.Subplot2D.plot(self.X[self.IndexSnap][k], self.Z[self.IndexSnap][k], marker='.', markersize=10, color=self.colorList[k])

            if self.indexRepres == 0:
                self.Subplot2D.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Z[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0)
            elif self.indexRepres == 1:
                self.Subplot2D.hist2d(self.X[self.IndexSnap][self.NbBodies-1:], self.Z[self.IndexSnap][self.NbBodies-1:], bins=[self.NbBinsX, self.NbBinsY], range=[[self.Xmin, self.Xmax],[self.Ymin, self.Ymax]], cmap='inferno')

            # Plot features
            self.Subplot2D.set_title('t='+str(round(self.t, 1))+' Myr')
            self.Subplot2D.set_xlabel('X (AU)')
            self.Subplot2D.set_xlim(self.Xmin, self.Xmax)
            self.Subplot2D.set_ylabel('Z (AU)')
            self.Subplot2D.set_ylim(self.Zmin, self.Zmax)
        
        elif self.indexView == 2: 

            self.RepresWidget.setEnabled(False)

            # Initialisation of axis
            self.Subplot3D = self.WindowPlot.Canvas.fig.add_subplot(111, projection='3d', aspect='equal')

            # Plot
            for k in range(self.NbBodies):
                # self.Subplot3D.plot(self.Xbf[k], self.Ybf[k], self.Zbf[k])
                self.Subplot3D.plot(self.Xbm[k], self.Ybm[k], self.Zbm[k], color=self.colorList[k], linestyle='--')
                if self.BodPosWidget.isChecked():
                    self.Subplot3D.plot(self.X[self.IndexSnap][k], self.Y[self.IndexSnap][k], self.Z[self.IndexSnap][k], marker='.', markersize=10, color=self.colorList[k])

            self.Subplot3D.scatter(self.X[self.IndexSnap][self.NbBodies-1:], self.Y[self.IndexSnap][self.NbBodies-1:], self.Z[self.IndexSnap][self.NbBodies-1:], s=self.SizePart, c='black', linewidths=0)

            # Plot features
            self.Subplot3D.set_title('t='+str(round(self.t, 1))+' Myr')
            self.Subplot3D.set_xlabel('X (AU)')
            self.Subplot3D.set_xlim(self.Xmin, self.Xmax)
            self.Subplot3D.set_ylabel('Y (AU)')
            self.Subplot3D.set_ylim(self.Ymin, self.Ymax)
            self.Subplot3D.set_zlabel('Z (AU)')
            self.Subplot3D.set_zlim(self.Zmin, self.Zmax)

        # Update canvas
        self.WindowPlot.Canvas.draw()



    def Ellipse2(self, a, e, Ex, Ey, Ez, Epx, Epy, Epz):

        E = linspace(-pi, pi, 100)
        x = Ex*a*(cos(E)-e) + Epx*a*sqrt(1-e**2)*sin(E)
        y = Ey*a*(cos(E)-e) + Epy*a*sqrt(1-e**2)*sin(E)
        z = Ez*a*(cos(E)-e) + Epz*a*sqrt(1-e**2)*sin(E)

        return x, y, z
    


class DiagramAE(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, e_m):
        super().__init__('Diagram a=f(e)', "Bodies' evolution of semi-major axis as a function of excentricity")

        # General data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.e_m = e_m

         # Plot initialisation
        self.Subplot = self.WindowPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        self.InitParams()


    def InitParams(self):
        
        # a limits
        self.Amin = 0
        self.AminWidget = LineEditInt('Amin', 'Semi-major axis minimum (AU)', self.Amin)
        self.WindowParam.Layout.addWidget(self.AminWidget)

        self.Amax = round(max(self.a_m[0]))
        self.AmaxWidget = LineEditInt('Amax', 'Semi-major axis maximum (AU)', self.Amax)
        self.AminWidget.Layout.addWidget(self.AmaxWidget) 

        # Particules' size
        self.SizePart = 0.05
        self.SizePartWidget = LineEditFloatPos('Size of particules', 'Plot size of particles', self.SizePart)
        self.WindowParam.Layout.addWidget(self.SizePartWidget)

        # Orbital resonance
        self.CheckRes = QCheckBox('Orbital resonance')
        self.WindowParam.Layout.addWidget(self.CheckRes)
        self.CheckRes.stateChanged.connect(self.ResStateChange)

        self.ButAddRes = QPushButton('+')
        self.WindowParam.Layout.addWidget(self.ButAddRes)
        self.ButAddRes.clicked.connect(self.AddResWidget)
        self.ButAddRes.setEnabled(False)

        self.PRes0 = 1

        self.PRef0 = 1

        self.ResWidgets = []
        self.ResWidgetsId = []
        self.NbWidgets0 = self.WindowParam.Layout.count()

        self.c = 0 # counter

    def AddResWidget(self):
        self.ResWidget = ResClass(self.NbBodies_m[self.IndexSnap], self.PRes0, self.PRef0)
        self.ResWidgets.append(self.ResWidget)
        self.c += 1
        self.ResWidget.Id = self.c
        self.ResWidgetsId.append(self.ResWidget.Id)
        self.WindowParam.Layout.addWidget(self.ResWidget)
        self.ResWidget.SignalDel.connect(self.DelResWidget)

    def DelResWidget(self, Id):
        index = self.ResWidgetsId.index(Id)
        self.ResWidgets.pop(index)
        self.ResWidgetsId.pop(index)

        for i in reversed(range(self.NbWidgets0, self.WindowParam.Layout.count())): 
            WidgetToRemove = self.WindowParam.Layout.itemAt(i).widget()
            self.WindowParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        
        for i in range(len(self.ResWidgets)):
            self.WindowParam.Layout.addWidget(self.ResWidgets[i])

        self.WindowParam.resize(self.WindowParam.minimumSize())


    def ResStateChange(self):
        ResState = self.CheckRes.isChecked()
        self.ButAddRes.setEnabled(ResState)
        for i in range(len(self.ResWidgets)):
            self.ResWidgets[i].setEnabled(ResState)




    # Update of parameters
    def UpdateParams(self):
        self.Amin = int(self.AminWidget.EditParam.text())
        self.Amax = int(self.AmaxWidget.EditParam.text())
        self.SizePart = float(self.SizePartWidget.EditParam.text().replace(',', '.'))
        if self.CheckRes.isChecked():
            self.nRefValues = []
            self.PResValues = []
            self.PRefValues = []
            for i in range(len(self.ResWidgets)):
                self.nRefValues.append(int(self.ResWidgets[i].nRefWidget.SpinParam.value()))
                self.PResValues.append(int(self.ResWidgets[i].PResWidget.EditParam.text()))
                self.PRefValues.append(int(self.ResWidgets[i].PRefWidget.EditParam.text()))
        

    # Plot
    def Plot(self):
        
        # Clear the plot
        self.Subplot.cla()

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.Canvas.draw()
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
        for k in range(self.NbBodies_m[self.IndexSnap]):
            self.Subplot.plot(self.a[k], self.e[k], color=self.colorList[k], marker='.', markersize=10)
        self.Subplot.scatter(self.a[self.NbBodies_m[self.IndexSnap]:], self.e[self.NbBodies_m[self.IndexSnap]:], s=self.SizePart, c='black', linewidths=0)

        if self.CheckRes.isChecked():
            for i in range(len(self.ResWidgets)):
                colorRef = self.colorList[self.nRefValues[i]]
                self.Subplot.axvline(self.aResValues[i], linewidth=1, color=colorRef, linestyle='--')
                # self.Subplot.annotate(str(self.PResValues[i])+':'+str(self.PRefValues[i]), xy=(self.aResValues[i], 0.9), bbox=dict(boxstyle='round', facecolor='white', edgecolor=colorRef))
                self.Subplot.text(self.aResValues[i], 0.9, ' '+str(self.PResValues[i])+':'+str(self.PRefValues[i])+' ', rotation=90, color=colorRef, bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=colorRef), fontsize=7, va='center', ha='center')

        # Plot features
        self.Subplot.set_title('t='+str(round(self.t, 1))+' Myr')
        self.Subplot.set_xlabel('Semi-major axis (AU)')
        self.Subplot.set_xlim(self.Amin, self.Amax)
        self.Subplot.set_ylabel('Eccentricity')
        self.Subplot.set_ylim(0, 1)
            
        # Update canvas
        self.WindowPlot.Canvas.draw()



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

        # Plot initialisation
        self.Subplot = self.WindowPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        self.InitParams()


    def InitParams(self):

        # Orbit number
        self.nOrbit = 0
        self.nOrbitWidget = SpinBox("Bodie's number",'Number of the bodie which is studied (counting from the center of the system outwards, including stars in last)', 0, self.NbBodies_f-1)
        self.WindowParam.Layout.addWidget(self.nOrbitWidget)

        # Orbit parameters
        self.ParamOrbitWidget = ComboBox('Orbital parameter', 'Orbit Parameter', ['a','e','i','W','w','M'])
        self.WindowParam.Layout.addWidget(self.ParamOrbitWidget)

        # Time limits
        self.Tmin = 0
        self.TminWidget = LineEditIntPos('Tmin', 'Time minimum (yr)', self.Tmin)
        self.WindowParam.Layout.addWidget(self.TminWidget)

        self.Tmax = round(max(self.t_f[0]))
        self.TmaxWidget = LineEditIntPos('Tmax', 'Time maximum (yr)', self.Tmax)
        self.TminWidget.Layout.addWidget(self.TmaxWidget) 

    def UpdateParams(self):
        self.nOrbit = int(self.nOrbitWidget.SpinParam.value())
        self.ParamOrbit = self.ParamOrbitWidget.ComboParam.currentText()
        self.Tmin = int(self.TminWidget.EditParam.text())
        self.Tmax = int(self.TmaxWidget.EditParam.text())


    # Plot
    def Plot(self):
        
        # Clear the plot
        self.Subplot.cla()

        # Update parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.Canvas.draw()
            return
        
        # Specific data
        t = [self.t_f[self.nOrbit],'Time (yr)']
        a = [self.a_f[self.nOrbit],'Semi-major axis (AU)']
        e = [self.e_f[self.nOrbit],'Eccentricity']
        i = [self.i[self.nOrbit],'Inclinaison (deg)']
        W = [self.W[self.nOrbit],'Longitude of ascending node (deg)']
        w = [self.w[self.nOrbit],'Argument of periastron (deg)']
        M = [self.M[self.nOrbit],'Initial mean longitude']
        
        # Plot with current parameters
        self.EvalParamOrbit = eval(self.ParamOrbit)
        self.Subplot.plot(t[0], self.EvalParamOrbit[0], linewidth=0.5)

        # Plot features
        self.Subplot.set_xlabel('Time (yr)')
        self.Subplot.set_xlim(self.Tmin, self.Tmax)
        self.Subplot.set_ylabel(self.EvalParamOrbit[1])
        
        # Update canvas
        self.WindowPlot.Canvas.draw()



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

        # Plot initialisation
        self.Subplot = self.WindowPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        self.InitParams()

    def InitParams(self):
        # Ordinate quantity formula
        self.YFormula = ''
        self.YFormulaWidget = LineEdit('Y formula','Ordinate quantity by combining t, a, e, i, w, W, [n] where n is the number of bodie (counting from the center of the system outwards, including stars in last), and math fonctions', self.YFormula)
        self.WindowParam.Layout.addWidget(self.YFormulaWidget)

        # Ordinate quantity label
        self.YLabel = ''
        self.YLabelWidget = LineEdit('Y label', 'Ordinate label', self.YLabel)
        self.WindowParam.Layout.addWidget(self.YLabelWidget)

        # Abscissa quantity formula
        self.XFormula = '' 
        self.XFormulaWidget = LineEdit('X formula','Abscissa quantity by combining t, a, e, i, w, W, [n] where n is the number of bodie (counting from the center of the system outwards, including stars in last), and math fonctions', self.XFormula)
        self.WindowParam.Layout.addWidget(self.XFormulaWidget)

        # Abscissa quantity label
        self.XLabel = ''
        self.XLabelWidget = LineEdit('X label', 'Abscissa label', self.XLabel)
        self.WindowParam.Layout.addWidget(self.XLabelWidget)

    def UpdateParameters(self):
        self.YFormula = self.YFormulaWidget.EditParam.text()
        self.YLabel = self.YLabelWidget.EditParam.text()
        self.XFormula = self.XFormulaWidget.EditParam.text()
        self.XLabel = self.XLabelWidget.EditParam.text()

    # Plot
    def Plot(self):

        # Clear the plot
        self.Subplot.cla()

        # Update of parameters
        try:
            self.UpdateParameters()
        except:
            print('Wrong Parameters')
            self.WindowPlot.Canvas.draw()
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
        
        # Update canvas
        self.WindowPlot.Canvas.draw()



class RadProfile(GeneralToolClass):
    def __init__(self, t_m, NbBodies_m, a_m, R):
        super().__init__('Radial profile', "Particules' integrated radial profile")

        # Data
        self.t_m = t_m
        self.NbBodies_m = NbBodies_m
        self.a_m = a_m
        self.R = R

        # Beta Pic other curves
        self.profileAug = transpose(loadtxt('/Users/lacquema/ByeGildas/Softs/SoftSwift/OtherCurves/bpic_Augereau_profile.dat', dtype = float))
        self.profileNE = transpose(loadtxt('/Users/lacquema/ByeGildas/Softs/SoftSwift/OtherCurves/bpic_Dent_profile_NE.dat', dtype = float))
        self.profileSW = transpose(loadtxt('/Users/lacquema/ByeGildas/Softs/SoftSwift/OtherCurves/bpic_Dent_profile_SW.dat', dtype = float))

        # Plot initialisation
        self.Subplot = self.WindowPlot.Canvas.fig.add_subplot(111)

        # Parameters initialisation
        self.InitParams()

    def InitParams(self):

        # R limits
        self.Rmin = 0
        self.RminWidget = LineEditIntPos('Rmin', 'Radii minimum (AU)', self.Rmin)
        self.WindowParam.Layout.addWidget(self.RminWidget)

        self.Rmax = round(max(self.R[0]))
        self.RmaxWidget = LineEditIntPos('Rmax', 'Radii maximum (AU)', self.Rmax)
        self.RminWidget.Layout.addWidget(self.RmaxWidget) 

        # Histogram
        self.Norm = ComboBox('Normalisation', 'Choice of the normalisation', ['None', 'One', 'Proportion'])
        self.Norm.ComboParam.setCurrentIndex(1)
        self.WindowParam.Layout.addWidget(self.Norm)

        self.NbBins = 100
        self.NbBinsWidget = LineEditIntPos('Bining', 'Number of bins', self.NbBins)
        self.WindowParam.Layout.addWidget(self.NbBinsWidget)

        # Check other curves
        self.CheckAugWidget = QCheckBox('Aug 2001')
        self.WindowParam.Layout.addWidget(self.CheckAugWidget)
        self.CheckAugWidget.setChecked(True)

        self.CheckDentWidget = QCheckBox('Dent 2014')
        self.WindowParam.Layout.addWidget(self.CheckDentWidget)
        self.CheckDentWidget.setChecked(True)

        self.CheckCurves = QCheckBox('Other curves')
        self.WindowParam.Layout.addWidget(self.CheckCurves)
        self.CheckCurves.stateChanged.connect(self.CurveStateChange)

        self.ButAddCurve = QPushButton('+')
        self.WindowParam.Layout.addWidget(self.ButAddCurve)
        self.ButAddCurve.clicked.connect(self.AddCurveWidget)
        self.ButAddCurve.setEnabled(False)

        self.CurveWidgets = []
        self.CurveWidgetsId = []
        self.NbWidgets0 = self.WindowParam.Layout.count()
        self.c = 0 # counter

    def AddCurveWidget(self):
        self.CurveWidget = CurveClass()
        self.CurveWidgets.append(self.CurveWidget)
        self.c += 1
        self.CurveWidget.Id = self.c
        self.CurveWidgetsId.append(self.CurveWidget.Id)
        self.WindowParam.Layout.addWidget(self.CurveWidget)
        self.CurveWidget.SignalDel.connect(self.DelCurveWidget)

    def DelCurveWidget(self, Id):
        index = self.CurveWidgetsId.index(Id)
        self.CurveWidgets.pop(index)
        self.CurveWidgetsId.pop(index)

        for i in reversed(range(self.NbWidgets0, self.WindowParam.Layout.count())): 
            WidgetToRemove = self.WindowParam.Layout.itemAt(i).widget()
            self.WindowParam.Layout.removeWidget(WidgetToRemove)
            WidgetToRemove.setParent(None)
        
        for i in range(len(self.CurveWidgets)):
            self.WindowParam.Layout.addWidget(self.CurveWidgets[i])

        self.WindowParam.resize(self.WindowParam.minimumSize())

    def CurveStateChange(self):
        CurveState = self.CheckCurves.isChecked()
        self.ButAddCurve.setEnabled(CurveState)
        for i in range(len(self.CurveWidgets)):
            self.CurveWidgets[i].setEnabled(CurveState)
    

    def UpdateParams(self):
        self.Rmin = int(self.RminWidget.EditParam.text())
        self.Rmax = int(self.RmaxWidget.EditParam.text())
        self.NbBins = int(self.NbBinsWidget.EditParam.text())
        if self.CheckCurves.isChecked():
            self.CurvePaths = []
            self.CurveLabels = []
            for i in range(len(self.CurveWidgets)):
                self.CurvePaths.append(self.CurveWidgets[i].PathWidget.EditParam.text())
                self.CurveLabels.append(self.CurveWidgets[i].LabelWidget.EditParam.text())
                if self.CurveLabels[i] == '': self.CurveLabels[i] = self.CurvePaths[i].split('/')[-1]
            

    # Plot
    def Plot(self):

        # Clear the plot
        self.Subplot.cla()

        # Update of parameters
        try:
            self.UpdateParams()
        except:
            print('Wrong Parameters')
            self.WindowPlot.Canvas.draw()
            return
        
        # Specific variables
        self.t = self.t_m[self.IndexSnap]/10**6

        # Plot with current parameters
        for k in range(self.NbBodies_m[self.IndexSnap]):
            self.Subplot.plot(self.a_m[self.IndexSnap][k], 0, color=self.colorList[k], marker='.', markersize=10)

        # Histogram
        histCount, histX = histogram([x for x in self.R[self.IndexSnap] if self.Rmin<x<self.Rmax], bins=self.NbBins)

        # Normalisation
        if self.Norm.ComboParam.currentIndex()==1:
            self.NormDiv = max(histCount)
        elif self.Norm.ComboParam.currentIndex()==2:
            self.NormDiv = sum(histCount)
        else:
            self.NormDiv = 1

        # Stairs
        self.Subplot.stairs(histCount/self.NormDiv, histX, label='Simulation', linewidth=1, color='black')
        
        # Other curves
        if self.CheckAugWidget.isChecked(): self.Subplot.plot(self.profileAug[0], self.profileAug[1], color='blue', linestyle='dashed', linewidth=0.5, label='Aug 2001')
        if self.CheckDentWidget.isChecked(): self.Subplot.plot((self.profileNE[0]+self.profileSW[0])/2, (self.profileNE[1]+self.profileSW[1])/2, color='red', linestyle='dashed', linewidth=0.5, label='Dent 2014')
        
        if self.CheckCurves.isChecked(): 
            try:
                for i in range(len(self.CurveWidgets)):
                    self.Curve = transpose(loadtxt(self.CurvePaths[i], dtype = float))
                    self.Subplot.plot(self.Curve[0], self.Curve[1], linestyle='dashed', linewidth=0.5, label=self.CurveLabels[i])
            except:
                print('Wrong file')
                print('File must be in the format:')
                print('Column 1 : abscissa')
                print('Column 2 : ordinate')

            

        # Plot features
        self.Subplot.legend()
        self.Subplot.set_title('t='+str(round(self.t, 1))+' Myr')
        self.Subplot.set_xlabel('Radius (AU)')
        self.Subplot.set_xlim(self.Rmin, self.Rmax)
        self.Subplot.set_ylabel('Number of particules')
        # self.Subplot.set_ylim(0, 1)
        
        # Update canvas
        self.WindowPlot.Canvas.draw()






### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    n=1
    PathFollowbodies = f'/Users/lacquema/Documents/Swiftdata/BetaPic/simu_bpicbcd_{n}/followbodies.dat'
    PathMextract = f'/Users/lacquema/Documents/Swiftdata/BetaPic/simu_bpicbcd_{n}/mextract.dat'       
    NbSteps, NbBodies_f, t_f, a_f, e_f, i, W, w, M = TransferDataClass.OpenFollowbodies(PathFollowbodies)
    NbSnapshots, t_m, NbBodies_m, NbParticles, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R = TransferDataClass.OpenMextract(PathMextract)
    ToolWidget = RadProfile(t_m, NbBodies_m, a_m, R)
    # ToolWidget = SpaceView(t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R)
    # ToolWidget = DiagramAE(t_m, NbBodies_m, a_m, e_m)
    ToolWidget.show() 
    app.exec() # Application execution
