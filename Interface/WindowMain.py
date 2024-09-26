#! /Users/lacquema/ByeGildas/bin/python3




##########################################################################################################
###             <<< Here you must input paths to mextract.dat and followbodies.dat >>>                 ###
##########################################################################################################
n=8
System = 'BetaPic'
PathFollowbodies = '/Users/lacquema/Documents/Swiftdata/'+System+f'/simu_bpicbcd_{n}/followbodies.dat'
PathMextract = '/Users/lacquema/Documents/Swiftdata/'+System+f'/simu_bpicbcd_{n}/mextract.dat'
##########################################################################################################




### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication, QProgressBar

# My packages
from SnapSelector import SnapSelectorClass
from Tools import *
from TransferData import TransferDataClass
from WindowMenu import LoadWindowClass



### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):

    def __init__(self, PathFollowbodies, PathMextract):
        super().__init__()

        # Data
        NbSteps, NbBodies_f, t_f, a_f, e_f, i, W, w, M = TransferDataClass.OpenFollowbodies(PathFollowbodies)
        NbSnapshots, t_m, NbBodies_m, NbParticles, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R = TransferDataClass.OpenMextract(PathMextract)

        # Window settings
        self.setWindowTitle('Swift data analysis of {}'.format(System))

        # Layout intialisation
        Layout = QVBoxLayout()

        # SnapSelector adding
        SnapSelectorWidget = SnapSelectorClass(NbSnapshots)
        Layout.addWidget(SnapSelectorWidget)

        # Separation
        Layout.addWidget(QProgressBar())

        # Space view tool adding
        SpaceViewWidget = SpaceView(t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R)
        Layout.addWidget(SpaceViewWidget)
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(SpaceViewWidget.Refresh_ActivePlots) # Refresh all active plots when the refresh button is clicked
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(SpaceViewWidget.Change_IndexSnap) # Refresh all active plots with new snapshot index when the value of this index is changed

        # Radial profile
        RadProfileWidget = RadProfile(t_m, NbBodies_m, a_m, R)
        Layout.addWidget(RadProfileWidget)
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(RadProfileWidget.Refresh_ActivePlots)
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(RadProfileWidget.Change_IndexSnap)

        # Diagram A=F(E)
        DiagramAEWidget = DiagramAE(t_m, NbBodies_m, a_m, e_m)
        Layout.addWidget(DiagramAEWidget)
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramAEWidget.Refresh_ActivePlots)
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(DiagramAEWidget.Change_IndexSnap)

        # Orbit Evolution tool adding
        DiagramTYWidget = DiagramTY(NbBodies_f, t_f, a_f, e_f, i, W, w, M)
        Layout.addWidget(DiagramTYWidget)
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramTYWidget.Refresh_ActivePlots)

        # Diagram Y=F(X)
        DiagramXYWidget = DiagramXY(NbBodies_f, t_f, a_f, e_f, i, W, w, M)
        Layout.addWidget(DiagramXYWidget)
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramXYWidget.Refresh_ActivePlots)

        # Widget container
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))
        
        # Showing
        self.show()
        LoadWin.close()

    # Close programme when the main window are closed
    def closeEvent(self, e):
        app.closeAllWindows()




if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    LoadWin = LoadWindowClass() # Loading window showing
    app.processEvents() # Continue the program
    WindowMain = WindowMainClass(PathFollowbodies, PathMextract) # Main window showing
    sys.exit(app.exec()) # Application execution
    