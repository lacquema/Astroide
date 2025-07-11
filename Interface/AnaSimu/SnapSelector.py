#! /Users/lacquema/Astroide.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtWidgets import QWidget, QHBoxLayout, QLabel, QSpinBox, QPushButton, QApplication
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import QTimer


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

        # Movie mode (play/pause) button
        self.BtnPlayPause = QPushButton()
        self.BtnPlayPause.setIcon(QIcon(f'{self.DirPath}/../Items/Userplay_button.png'))
        self.BtnPlayPause.setStatusTip('Play/Stop movie')
        self.isPlaying = False
        Layout.addWidget(self.BtnPlayPause)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.Auto_next_snap)

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

        # Input number for maximum play time
        self.EditTimePlay = QSpinBox() 
        self.EditTimePlay.setRange(1, 120) # Input number for maximum play time (in seconds)
        self.EditTimePlay.setValue(10)  # default value
        self.EditTimePlay.setFixedWidth(50)  
        Layout.addWidget(self.EditTimePlay)
        self.Update_play_interval() # Set play interval based on EditTimePlay value
        self.EditTimePlay.valueChanged.connect(self.Update_play_interval)
        self.EditTimePlay.setStatusTip('Play time [s]')

        # Click actions
        self.BtnFirstSnap.clicked.connect(self.GoFirstSnap)
        self.BtnPreviousSnap.clicked.connect(self.GoPreviousSnap) # Ask to previous snapshot
        self.BtnPlayPause.clicked.connect(self.Toggle_play_stop) # Ask to play or pause
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

    def Update_play_interval(self):
        # Total play time in seconds from EditTimePlay
        total_play_time_s = self.EditTimePlay.value()
        # Calculate interval in ms between snapshots
        if self.NbSnap > 0:
            self.play_interval_ms = int((total_play_time_s * 1000) / self.NbSnap)
        else:
            self.play_interval_ms = 1000
        # If playing, update timer interval immediately
        if hasattr(self, 'timer') and self.isPlaying:
            self.timer.stop()
            self.timer.start(self.play_interval_ms)

    # Toggle play/pause for movie mode
    def Toggle_play_stop(self):
        if not self.isPlaying:
            # Start playing
            self.isPlaying = True
            self.BtnPlayPause.setIcon(QIcon(f'{self.DirPath}/../Items/Pause.png'))
            # If at last snapshot, restart from 0
            if self.IndexSnap == self.NbSnap:
                self.EditIndexSnap.setValue(0)
            self.timer.start(self.play_interval_ms)
        else:
            # Pause
            self.isPlaying = False
            self.BtnPlayPause.setIcon(QIcon(f'{self.DirPath}/../Items/Userplay_button.png'))
            self.BtnPlayPause.setText('')
            self.timer.stop()

    # Advance to next snapshot automatically
    def Auto_next_snap(self):
        if self.IndexSnap < self.NbSnap:
            self.EditIndexSnap.setValue(self.IndexSnap + 1)
        else:
            self.Toggle_play_stop()  # Stop at the end

# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    SnapSelectorWidget = SnapSelectorClass(20, 20000000)
    SnapSelectorWidget.show()
    app.exec() # Application execution
