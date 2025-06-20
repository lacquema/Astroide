#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QPushButton, QCheckBox
from PyQt6.QtCore import pyqtSignal

# My packages
from Parameters import *


class ResClass(QWidget):

    SignalDel = pyqtSignal(int) # Signal to transfer the identity of the resonance

    def __init__(self, NbBodies):
        super().__init__()

        # Parameters
        self.NbBodies = NbBodies  # Number of bodies in the system

        # Identity
        self.Id = 0
        
        # Layout
        self.Layout = QHBoxLayout()

        # Reference number of the resonance
        self.ListBody = [str(k) for k in range(self.NbBodies)][1:]
        self.nRefWidget = ComboBox(None, 'Number of the orbit which is the reference of the resonance counted in the order of the input file', self.ListBody)
        self.Layout.addWidget(self.nRefWidget)

        # Period ratio between the two resonating bodies
        self.PResWidget = SpinBox(None, 'Resonant orbit period', 1, 1, None)
        self.Layout.addWidget(self.PResWidget)

        self.LabelRatioWidget = QLabel(' :')
        self.Layout.addWidget(self.LabelRatioWidget)

        self.PRefWidget = SpinBox(None, 'Reference orbit period', 1, 1, None)
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