#! /Users/lacquema/Astroide.env/bin/python3
import sys
import os

### --- Packages --- ###

# Transverse packages

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication, QVBoxLayout, QPushButton, QFileDialog
from PyQt6.QtCore import Qt, pyqtSignal

# My packages
# from NewSimu.TabsParamNew import *
from Parameters import *
from UtilsAnaSimu import DelAllWidgetsBtw
from WindowMain import WindowMainClass
# from WindowMenu import LoadWindowClass


### --- Parameters Window Generating --- ###

class WindowSetAnaSimu(QMainWindow):

    SignalCloseWindowSetAnaSimu = pyqtSignal()
    ReSignalCloseWindowMain = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Window characteristics
        self.setWindowTitle('Settings of the analysis')
        self.setMinimumWidth(600)

        # Layout initialisation
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset all tab settings')
        self.Layout.addWidget(self.BtnReset)
        self.BtnReset.clicked.connect(self.ResetParams)
        
        self.InitWidgets()

        # Widget Container
        self.Container = QWidget()
        self.Container.setLayout(self.Layout)

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))



    def InitWidgets(self):
        
        self.SimuPath = PathBrowser('Directory path', 'Path to the data to analyse', 0)
        self.Layout.addWidget(self.SimuPath)

        self.FollowbodiesFileName = LineEdit('Followbodies file', 'Name of the followbodies file with extension', 'followbodies.dat')
        self.Layout.addWidget(self.FollowbodiesFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.MextractFileName = LineEdit('Mextract file', 'Name of the mextract file with extension', 'mextract.dat')
        self.Layout.addWidget(self.MextractFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()

    def AnalyseSimu(self):
        if len(self.SimuPath.EditPath.text()) == 0:
            print('Simulation directory path not given.')
            print('Check your inputs.')
        else:
            self.OpenWinMain()


    def OpenWinMain(self):
        try:
            self.BtnStart.setEnabled(False)
            self.WinMain = WindowMainClass(self.SimuPath.EditPath.text()+self.FollowbodiesFileName.EditParam.text(), self.SimuPath.EditPath.text()+self.MextractFileName.EditParam.text())
            self.WinMain.SignalCloseWindowMain.connect(self.ReSignalCloseWindowMain.emit)
            self.WinMain.show()
            self.close()
        except:
            print('Data not found: check the directory path and the name of input files.')

    # def FindInputFiles(self):
    #     Files = os.listdir(self.SimuPath.EditPath.text()[:-1])
    #     self.SimuName.EditParam.setText(self.SimuPath.EditPath.text()[:-1].split('/')[-1])
    #     txtFiles = []
    #     datFiles = []
    #     for x in Files:
    #         if x[-4:] == '.txt':
    #             txtFiles.append(x)
    #         if x[-4:] == '.dat':
    #             datFiles.append(x)
    #     if 'logfile.dat' in datFiles: datFiles.remove('logfile.dat')
    #     if len(txtFiles) == 1:
    #         self.DataFileName.EditParam.setText(txtFiles[0])
    #     else:
    #         self.DataFileName.EditParam.setText('')
    #         print('Data file not found')
    #     if len(datFiles) == 1:
    #         self.SimuFileName.EditParam.setText(datFiles[0])
    #     else:
    #         self.SimuFileName.EditParam.setText('')
    #         print('Simulation file not found')
        

    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        try: # In the case self.WinMain openning, dont show again WindowMenu
            if self.WinMain.isVisible() == False:
                self.SignalCloseWindowSetAnaSimu.emit() 
        except:
            self.SignalCloseWindowSetAnaSimu.emit() 
        




        


# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetAnaSimu()
    WindowParam.show()
    app.exec() # Application execution
