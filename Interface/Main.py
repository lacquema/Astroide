#! /Users/lacquema/Astroide.env/bin/python3

import sys
import os

# Adding necessary paths to import custom modules
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.dirname(__file__) + '/AnaSimu')
sys.path.append(os.path.dirname(__file__) + '/NewSimu')
sys.path.append(os.path.dirname(__file__) + '/ContSimu')

### --- Importing packages --- ###

# PyQt packages
from PyQt6.QtWidgets import QApplication

# Importing custom window classes
from WindowMenu import WindowMenuClass
from NewSimu.WindowSetNewSimu import WindowSetNewSimu
from ContSimu.WindowSetContSimu import WindowSetContSimu
from AnaSimu.WindowSetAnaSimu import WindowSetAnaSimu


### --- Main class to manage windows --- ###

class MainClass:

    def __init__(self):
        super().__init__()

        # Initialize the windows
        self.WinMenu = WindowMenuClass()  # Main menu window
        self.WinSetNewSimu = WindowSetNewSimu()  # Window for starting a new simulation
        self.WinSetContSimu = WindowSetContSimu()  # Window for continuing a simulation
        self.WinSetAnaSimu = WindowSetAnaSimu()  # Window for analyzing a simulation

        # Connect signals to handle window closures
        self.WinSetNewSimu.SignalCloseWindowSetNewSimu.connect(self.ReOpenWinLoad)
        self.WinSetContSimu.SignalCloseWindowSetContSimu.connect(self.ReOpenWinLoad)
        self.WinSetAnaSimu.SignalCloseWindowSetAnaSimu.connect(self.ReOpenWinLoad)
        self.WinSetAnaSimu.ReSignalCloseWindowMain.connect(self.ReOpenWinLoad)

        # Show the main menu window
        self.WinMenu.show()

        # Connect main menu buttons to their respective actions
        self.WinMenu.BtnNew.clicked.connect(self.OpenWinSetNewSimu)  # Open new simulation window
        self.WinMenu.BtnContinue.clicked.connect(self.OpenWinSetContSimu)  # Open continue simulation window
        self.WinMenu.BtnAnalyse.clicked.connect(self.OpenWinSetAnaSimu)  # Open analyze simulation window

    def OpenWinSetNewSimu(self):
        """
        Method to open the new simulation window
        """
        self.WinMenu.close()  # Close the main menu
        self.WinSetNewSimu.show()  # Show the new simulation window

    def OpenWinSetContSimu(self):
        """
        Method to open the continue simulation window
        """
        self.WinMenu.close()  # Close the main menu
        self.WinSetContSimu.show()  # Show the continue simulation window

    def OpenWinSetAnaSimu(self):
        """
        Method to open the analyze simulation window
        """
        self.WinMenu.close()  # Close the main menu
        self.WinSetAnaSimu.show()  # Show the analyze simulation window

    def ReOpenWinLoad(self):
        """
        Method to reopen the main menu window
        """
        app.closeAllWindows()  # Close all open windows
        self.WinMenu.show()  # Show the main menu again

    def closeEvent(self, e):
        """
        Method to handle the closure of the main window
        """
        self.WinMenu.show()  # Show the main menu again if the window is closed


### --- Main entry point of the application --- ###

if __name__ == "__main__":
    app = QApplication(sys.argv)  # Create the PyQt application
    WindowMain = MainClass()  # Instantiate the main class
    sys.exit(app.exec())  # Execute the application

