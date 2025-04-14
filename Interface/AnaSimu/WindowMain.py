#! /Users/lacquema/Astroide.env/bin/python3

### --- Packages --- ###

# Standard Python packages
import sys
import os
import numpy as np
from scipy.interpolate import interp1d

# PyQt packages for GUI components and signals
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QStatusBar, QWidget, QApplication, QProgressBar, QGridLayout
from PyQt6.QtCore import pyqtSignal

# Custom modules
from SnapSelector import SnapSelectorClass
from Tools import *
from TransferData import TransferDataClass
from WindowLoad import LoadWindowClass


### --- Main Window Generating --- ###

class WindowMainClass(QMainWindow):
    # Signal emitted when the main window is closed
    SignalCloseWindowMain = pyqtSignal()


    def __init__(self, PathFollowbodies=str, PathMextract=str):
        super().__init__()

        # Load data from the provided file paths
        NbSteps, NbBodies_f, t_f, a_f, e_f, i, W, w, M = TransferDataClass.OpenFollowbodies(PathFollowbodies)
        NbSnapshots, t_m, NbBodies_m, NbParticles, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R = TransferDataClass.OpenMextract(PathMextract)
        tmax = np.max(t_m)  # Maximum time value for the simulation

        # Set the window title based on the simulation folder name
        self.setWindowTitle('Astroide analysis of {}'.format(PathFollowbodies.split('/')[-2]))

        # Initialize the main layout as a vertical box layout
        Layout = QVBoxLayout()

        # Add the SnapSelector widget for snapshot selection
        SnapSelectorWidget = SnapSelectorClass(NbSnapshots-1, tmax)
        Layout.addWidget(SnapSelectorWidget)

        # Add a visual delimiter for separation
        Layout.addWidget(Delimiter(Title='Plots :'))

        # Grid layout
        GridLayout = QGridLayout()

        # Add the SpaceView widget for 3D space visualization
        SpaceViewWidget = SpaceView(t_m, NbBodies_m, a_m, e_m, Ex, Ey, Ez, Epx, Epy, Epz, X, Y, Z, R)
        GridLayout.addWidget(SpaceViewWidget, 0, 0, 1, 3)

        # Connect SnapSelector signals to SpaceView methods for interactivity
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(SpaceViewWidget.Refresh_active_plots)  # Refresh plots on button click
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(SpaceViewWidget.Change_IndexSnap)  # Update plots on index change

        # Add the Radial Profile widget for radial distribution visualization
        RadProfileWidget = RadProfile(t_m, NbBodies_m, a_m, X, Y, Z, R)
        GridLayout.addWidget(RadProfileWidget, 0, 3, 1, 3)

        # Connect SnapSelector signals to RadProfile methods for interactivity
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(RadProfileWidget.Refresh_active_plots)
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(RadProfileWidget.Change_IndexSnap)

        # Add the Diagram A=F(E) widget for energy vs semi-major axis visualization
        DiagramAEWidget = DiagramAE(t_m, NbBodies_m, a_m, e_m)
        GridLayout.addWidget(DiagramAEWidget, 2, 0, 1, 2)

        # Connect SnapSelector signals to DiagramAE methods for interactivity
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramAEWidget.Refresh_active_plots)
        SnapSelectorWidget.EditIndexSnap.valueChanged.connect(DiagramAEWidget.Change_IndexSnap)

        # Add the Orbit Evolution widget for orbital parameter evolution visualization
        DiagramTYWidget = DiagramTY(NbBodies_f, t_f, a_f, e_f, i, W, w, M)
        GridLayout.addWidget(DiagramTYWidget, 2, 2, 1, 2)

        # Connect SnapSelector signals to DiagramTY methods for interactivity
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramTYWidget.Refresh_active_plots)

        # Add the Diagram Y=F(X) widget for position-based visualization
        DiagramXYWidget = DiagramXY(NbBodies_f, t_f, a_f, e_f, i, W, w, M)
        GridLayout.addWidget(DiagramXYWidget, 2, 4, 1, 2)

        # Connect SnapSelector signals to DiagramXY methods for interactivity
        SnapSelectorWidget.BtnRefreshSnap.clicked.connect(DiagramXYWidget.Refresh_active_plots)

        # Add grid layout to layout
        Layout.addLayout(GridLayout)

        # Create a container widget to hold the layout and set it as the central widget
        Container = QWidget()
        Container.setLayout(Layout)
        self.setCentralWidget(Container)

        # Add a status bar to the main window
        self.setStatusBar(QStatusBar(self))

        # Display the main window
        self.show()

    # Override the closeEvent to handle application closure
    def closeEvent(self, e):
        try:
            # Attempt to close all application windows
            app.closeAllWindows()
        except:
            # Emit the custom signal if an exception occurs
            self.SignalCloseWindowMain.emit()


# Entry point of the application
if __name__ == "__main__":
    app = QApplication(sys.argv)  # Create the application instance
    # LoadWin = LoadWindowClass()  # Uncomment to show a loading window
    # app.processEvents()  # Process events to continue the program
    PathFollowbodies = '/Users/lacquema/Simulations/twa7/twa7_cb_a_dyn/twa7_cb_a_dyn_2/followbodies.dat'
    PathMextract = '/Users/lacquema/Simulations/twa7/twa7_cb_a_dyn/twa7_cb_a_dyn_2/mextract.dat'
    WindowMain = WindowMainClass(PathFollowbodies, PathMextract)  # Create and show the main window
    sys.exit(app.exec())  # Execute the application
