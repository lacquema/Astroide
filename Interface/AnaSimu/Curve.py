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

        # Browse button
        self.BtnBrowse = QPushButton('Browse')
        self.Layout.addWidget(self.BtnBrowse)
        self.BtnBrowse.clicked.connect(self.DialBrowse)

        # Path of the curve file
        self.PathFile = ''

        # Name of the curve file
        self.NameFile = LineEdit(None, 'Name of the added curve file', '')
        self.NameFile.EditParam.setMinimumWidth(200)
        # self.NameFile.setEnabled(False)
        self.NameFile.EditParam.setReadOnly(True)
        self.Layout.addWidget(self.NameFile)

        # # Label of the curve
        # self.LabelWidget = LineEdit('Label', 'Label of the added curve', '')
        # self.Layout.addWidget(self.LabelWidget)

        # Delete the curve
        self.ButtonDel = QPushButton('-')
        self.Layout.addWidget(self.ButtonDel)
        self.ButtonDel.clicked.connect(lambda: self.SignalDel.emit(self.Id))

        # Widget container  
        self.setLayout(self.Layout)


    def DialBrowse(self):
        self.File = QFileDialog.getOpenFileName(self)
        self.PathFile = self.File[0]
        self.NameFile.EditParam.setText(self.PathFile.split('/')[-1])


### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    ResWidget = CurveClass()
    ResWidget.show()
    app.exec() # Application execution