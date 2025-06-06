#! /Users/lacquema/Astroide.env/bin/python3
import sys
import os
sys.path.append(os.path.dirname(__file__)+'/..')

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
from WindowWithFinder import WindowWithFinder
# from WindowMenu import LoadWindowClass


### --- Parameters Window Generating --- ###

class WindowSetAnaSimu(WindowWithFinder):

    SignalCloseWindowSetAnaSimu = pyqtSignal()
    ReSignalCloseWindowMain = pyqtSignal()
    
    def __init__(self):
        super().__init__()

        # Window characteristics
        self.setWindowTitle('Settings of the analysis')
        self.setMinimumWidth(1000)
        self.setMinimumHeight(500)

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

        # Add container to the split main window
        self.Splitter.addWidget(self.Container)

        # Connect folder to edit path
        self.Finder.doubleClicked.connect(self.ChangePath)

        # Status bar
        self.setStatusBar(QStatusBar(self))


    def InitWidgets(self):
        
        self.SimuFilePathW = LineEdit('Path', 'Path to the simulation files to analyse', '')
        self.Layout.addWidget(self.SimuFilePathW)

        self.FollowbodiesFileName = LineEdit('Followbodies file', 'Name of the followbodies file with extension', '')
        self.Layout.addWidget(self.FollowbodiesFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.MextractFileName = LineEdit('Mextract file', 'Name of the mextract file with extension', '')
        self.Layout.addWidget(self.MextractFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.BtnStart = QPushButton('Analyse the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)
        self.BtnStart.clicked.connect(self.AnalyseSimu)

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)

    def ChangePath(self):
        index = self.Finder.selectedIndexes()[0]
        info = self.Model.fileInfo(index)
        file_path = info.absoluteFilePath()
        if self.is_valid_file(file_path):  # Check if the file is valid
            self.SimuFilePathW.EditParam.setText(file_path + '/')
            # Set Followbodies file if exists
            followbodies_path = os.path.join(file_path, 'followbodies.dat')
            if os.path.isfile(followbodies_path):
                self.FollowbodiesFileName.EditParam.setText('followbodies.dat')
            else:
                self.FollowbodiesFileName.EditParam.setText('')
                print('\nFollowbodies file not found, it may has been renamed. Please check the directory.')
            # Set Mextract file if exists
            mextract_path = os.path.join(file_path, 'mextract.dat')
            if os.path.isfile(mextract_path):
                self.MextractFileName.EditParam.setText('mextract.dat')
            else:
                self.MextractFileName.EditParam.setText('')
                print('\nMextract file not found, it may has been renamed. Please check the directory.')
        else:
            print('\nSelected directory is not valid')
            self.ClearEdits()

    def is_valid_file(self, file_path):
        # Check if the given path is a directory
        if os.path.isdir(file_path):
            # Check if the directory contains at least one .dat file
            for file_name in os.listdir(file_path):
                if file_name.endswith('.dat'):
                    return True
        return False

    def ClearEdits(self):
        # self.DataFileNameW.EditParam.setText('')
        self.SimuFilePathW.EditParam.setText('')
        self.FollowbodiesFileName.EditParam.setText('')
        self.MextractFileName.EditParam.setText('')
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()

    def AnalyseSimu(self):
        if len(self.SimuFilePathW.EditParam.text()) == 0:
            print('\nSimulation path not given')
        else:
            try:
                self.OpenWinMain()
            except Exception as e:
                print(f'\nThere is a problem with inputs: {e}')

    def OpenWinMain(self):
        self.WinMain = WindowMainClass(self.SimuFilePathW.EditParam.text()+self.FollowbodiesFileName.EditParam.text(), self.SimuFilePathW.EditParam.text()+self.MextractFileName.EditParam.text())
        self.WinMain.SignalCloseWindowMain.connect(self.ReSignalCloseWindowMain.emit)
        self.close()
        self.WinMain.show()

        

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
