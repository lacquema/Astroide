#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QPushButton, QWidget, QStatusBar, QApplication


### --- Parameters Window Generating --- ###

class WindowParamClass(QMainWindow):
    
    SignalCloseWindowParam = pyqtSignal() # initiation of the closeEvent signal
    
    def __init__(self, ToolName):
        super().__init__()

        # Window characteristics
        self.setWindowTitle(ToolName+': Parameters')

        # Layout
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset')
        self.Layout.addWidget(self.BtnReset)

        # Widget container
        self.Container = QWidget()
        self.Container.setLayout(self.Layout)
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowParam.emit() 



# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowParamClass('ToolName')
    WindowParam.resize(250, 200)
    WindowParam.show()
    app.exec() # Application execution