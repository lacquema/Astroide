#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QPushButton, QCheckBox, QFileDialog
from PyQt6.QtCore import pyqtSignal

# My packages
from Parameters import *


class CurveClass(QWidget):

    SignalDel = pyqtSignal(int) # Signal to transfer the identity of the resonance

    def __init__(self):
        super().__init__()

        # Identity
        self.Id = 0
        
        # Layout
        self.Layout = QHBoxLayout()

        # Path of the curve
        self.PathWidget = LineEdit('Path', 'Path of the added curve', '')
        self.Layout.addWidget(self.PathWidget)

        # Browse button
        self.BtnBrowse = QPushButton('Browse')
        self.Layout.addWidget(self.BtnBrowse)
        self.BtnBrowse.clicked.connect(self.DialBrowse)

        # Label of the curve
        self.LabelWidget = LineEdit('Label', 'Label of the added curve', '')
        self.Layout.addWidget(self.LabelWidget)

        # Delete the curve
        self.ButtonDel = QPushButton('-')
        self.Layout.addWidget(self.ButtonDel)
        self.ButtonDel.clicked.connect(lambda: self.SignalDel.emit(self.Id))

        # Widget container  
        self.setLayout(self.Layout)


    def DialBrowse(self):
        self.File = QFileDialog.getOpenFileName(self)
        self.PathWidget.EditParam.setText(self.File[0])


### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    ResWidget = CurveClass()
    ResWidget.show()
    app.exec() # Application execution