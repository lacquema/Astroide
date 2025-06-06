#! /Users/lacquema/Astroide.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QSpinBox, QPushButton, QApplication
from PyQt6.QtGui import QIcon


### --- Snap Selector Widget Genereting --- ###

class SnapSelectorClass(QWidget):
    
    def __init__(self, NbSnap, tmax):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Variables
        self.NbSnap = NbSnap
        # print(self.NbSnap)
        self.IndexSnap = 0
        self.tmax = tmax
        
        # Layout
        Layout = QHBoxLayout()

        # Label snapshot number
        LblIndexSnap = QLabel("Snapshot number :")
        Layout.addWidget(LblIndexSnap)

        # Input number
        self.EditIndexSnap = QSpinBox()
        self.EditIndexSnap.setRange(0, NbSnap)
        Layout.addWidget(self.EditIndexSnap)

        Layout.addSpacing(5)

        # Label time
        self.LblTime = QLabel()
        self.LblTime.setFixedWidth(120)
        self.ChangeTime()
        Layout.addWidget(self.LblTime)

        # Button to fist snapshot
        self.BtnFirstSnap = QPushButton()
        self.BtnFirstSnap.setIcon(QIcon(f'{self.DirPath}/../Items/arrowStopL.png'))
        self.BtnFirstSnap.setEnabled(False)
        self.BtnFirstSnap.setStatusTip('First snapshot')
        Layout.addWidget(self.BtnFirstSnap)
        
        # Button to previous snapshot
        self.BtnPreviousSnap = QPushButton()
        self.BtnPreviousSnap.setIcon(QIcon(f'{self.DirPath}/../Items/arrowL.png'))
        self.BtnPreviousSnap.setEnabled(False)
        self.BtnPreviousSnap.setStatusTip('Previous snapshot')
        Layout.addWidget(self.BtnPreviousSnap)
        
        # Button to refresh snapshot
        self.BtnRefreshSnap = QPushButton()
        self.BtnRefreshSnap.setIcon(QIcon(f'{self.DirPath}/../Items/arrowCircle.png'))
        self.BtnRefreshSnap.setStatusTip('Refresh')
        Layout.addWidget(self.BtnRefreshSnap)

        # Button to next snapshot
        self.BtnNextSnap = QPushButton()
        self.BtnNextSnap.setIcon(QIcon(f'{self.DirPath}/../Items/arrowR.png'))
        self.BtnNextSnap.setStatusTip('Next snapshot')
        Layout.addWidget(self.BtnNextSnap)

        # Button to last snapshot
        self.BtnLastSnap = QPushButton()
        self.BtnLastSnap.setIcon(QIcon(f'{self.DirPath}/../Items/arrowStopR.png'))
        self.BtnLastSnap.setStatusTip('Last snapshot')
        Layout.addWidget(self.BtnLastSnap)

        # Click actions
        self.BtnFirstSnap.clicked.connect(self.GoFirstSnap)
        self.BtnPreviousSnap.clicked.connect(self.GoPreviousSnap) # Ask to previous snapshot
        self.BtnNextSnap.clicked.connect(self.GoNextSnap) # Ask to next snapshot
        self.BtnLastSnap.clicked.connect(self.GoLastSnap)
        self.EditIndexSnap.valueChanged.connect(self.ChangeValueIndexSnap) # Enable button on limits

        # Widget container
        self.setLayout(Layout) # SnapSelector is directly the widget container

    # Change index snapshot and state of button
    def ChangeValueIndexSnap(self, value):
        self.IndexSnap = value
        self.ChangeTime()
        if self.IndexSnap == 0:
            self.BtnFirstSnap.setEnabled(False)
            self.BtnPreviousSnap.setEnabled(False)
            self.BtnNextSnap.setEnabled(True)
            self.BtnLastSnap.setEnabled(True)
        elif self.IndexSnap == self.NbSnap:
            self.BtnFirstSnap.setEnabled(True)
            self.BtnPreviousSnap.setEnabled(True)
            self.BtnNextSnap.setEnabled(False)
            self.BtnLastSnap.setEnabled(False)
        else :
            self.BtnFirstSnap.setEnabled(True)
            self.BtnPreviousSnap.setEnabled(True)
            self.BtnNextSnap.setEnabled(True)
            self.BtnLastSnap.setEnabled(True)

    # Change time label
    def ChangeTime(self):
        self.time = self.tmax/self.NbSnap*self.IndexSnap/10**6
        self.LblTime.setText(f"=>   t = {self.time:.1f} Myr")

    # Change index snapshot to the first directly
    def GoFirstSnap(self):
        self.IndexSnap = 0
        self.EditIndexSnap.setValue(self.IndexSnap)

    # Change index snapshot to previous
    def GoPreviousSnap(self):
        if self.IndexSnap != 0:
            self.IndexSnap -= 1
            self.EditIndexSnap.setValue(self.IndexSnap)

    # Change index snapshot to next
    def GoNextSnap(self):
        if self.IndexSnap != self.NbSnap:
            self.IndexSnap += 1
            self.EditIndexSnap.setValue(self.IndexSnap)

    # Change index snapshot to the last directly
    def GoLastSnap(self):
        self.IndexSnap = self.NbSnap
        self.EditIndexSnap.setValue(self.IndexSnap)

# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    SnapSelectorWidget = SnapSelectorClass(20, 20000000)
    SnapSelectorWidget.show()
    app.exec() # Application execution
