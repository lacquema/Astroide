#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QPushButton, QCheckBox
from PyQt6.QtCore import pyqtSignal

# My packages
from Parameters import *


class ResClass(QWidget):

    SignalDel = pyqtSignal(int) # Signal to transfer the identity of the resonance

    def __init__(self, NbBodies, PRes, PRef):
        super().__init__()

        # Identity
        self.Id = 0
        
        # Layout
        self.Layout = QHBoxLayout()

        # Reference number of the resonance
        self.nRefWidget = SpinBox(None, 'Number of the bodie which is the reference of the resonance (counting from the center of the system outwards, including stars in first)', 0, NbBodies-1)
        self.Layout.addWidget(self.nRefWidget)

        # Period ratio between the two resonating bodies
        self.PResWidget = LineEditIntPos(None, 'Resonant orbit period', PRes)
        self.Layout.addWidget(self.PResWidget)

        self.LabelRatioWidget = QLabel(' /')
        self.Layout.addWidget(self.LabelRatioWidget)

        self.PRefWidget = LineEditIntPos(None, 'Reference orbit period', PRef)
        self.Layout.addWidget(self.PRefWidget)

        # Delete the resonance
        self.ButtonDel = QPushButton('-')
        self.Layout.addWidget(self.ButtonDel)
        self.ButtonDel.clicked.connect(lambda: self.SignalDel.emit(self.Id))

        # Widget container
        self.setLayout(self.Layout)

    


### --- Check --- ###
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    ResWidget = ResClass(3, 1, 1)
    ResWidget.show()
    app.exec() # Application execution