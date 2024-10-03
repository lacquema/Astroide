#! /Users/lacquema/Oracle.env/bin/python3
import sys
import os
import shutil
import subprocess
sys.path.append(os.path.dirname(__file__))


### --- Packages --- ###

# Transverse packages
import numpy as np

# PyQt packages
from PyQt6.QtWidgets import QTabWidget, QMainWindow, QStatusBar, QApplication

# My packages
from Tabs import *


### --- Parameters Window Generating --- ###

class WindowSetNewSimu(QMainWindow):

    SignalCloseWindowSetNewSimu = pyqtSignal() # initiation of the closeEvent signal
    
    def __init__(self):
        super().__init__()

        # Directory path
        self.DirPath = os.path.dirname(__file__)

        # Window characteristics
        self.setWindowTitle('Settings of the new simulation')

        # Widget Container
        self.Container = QTabWidget()

        # Tab 1
        self.TabSimuFiles = TabSimuFiles()
        self.Container.addTab(self.TabSimuFiles, 'Simulation files')
        self.TabSimuFiles.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)
        self.TabSimuFiles.SimuName.EditParam.textChanged.connect(self.ChangeStartOrder)

        # Tab 2
        self.TabSimuSets = TabSimuSets()
        self.Container.addTab(self.TabSimuSets, 'Simulation settings')

        # Tab 3
        self.TabOrbitsParams = TabOrbitsParams()
        self.Container.addTab(self.TabOrbitsParams, 'Orbital parameters')
        self.TabOrbitsParams.NbPart.SpinParam.valueChanged.connect(self.ChangeSumaPara)

        # Tab 4
        self.TabStart = TabStart()
        self.Container.addTab(self.TabStart, 'Starting')
        self.TabStart.NbHours.SpinParam.valueChanged.connect(self.ChangeStartOrder)
        self.TabStart.NbCores.SpinParam.valueChanged.connect(self.ChangeStartOrder)
        self.TabStart.NbCores.SpinParam.valueChanged.connect(self.ChangeSumaPara)
        self.TabStart.NbBranch.SpinParam.valueChanged.connect(self.ChangeSumaPara)
        self.ChangeSumaPara()
        self.TabStart.BtnReset.clicked.connect(self.ChangeSumaPara)
        self.TabStart.BtnCreate.clicked.connect(self.CreateInputFiles)
        # self.TabStart.BtnStart.clicked.connect(self.StartSimu)
        

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))


    def ChangeSumaPara(self):
        self.NbCoresPerSimu = self.TabStart.NbCores.SpinParam.value()//self.TabStart.NbBranch.SpinParam.value()
        self.NbPartPerSimu = int(np.ceil(self.TabOrbitsParams.NbPart.SpinParam.value()/self.TabStart.NbBranch.SpinParam.value()))
        self.TabStart.SumaPara.setText(f'=> {self.NbPartPerSimu} particules and {self.NbCoresPerSimu} cores per simulation')


    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetNewSimu.emit() 

    def ChangeStartOrder(self):
        self.TabStart.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core={self.NbCoresPerSimu},walltime={self.TabStart.NbHours.SpinParam.value()} --project dynapla {self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()}/{self.TabSimuFiles.InGenFileName.EditParam.text()}')
    
    def CreateInputFiles(self):
        print('--------------------------')
        if len(self.TabSimuFiles.SimuPath.EditPath.text()) == 0 or len(self.TabSimuFiles.SimuName.EditParam.text()) == 0:
            print('Simulation path not given.')
            print('Check your inputs.')
        else:
            if os.path.exists(self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()):
                print('This directory already exists.')
                print('Check your inputs.')
            else: 
                os.makedirs(self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text())
                print(f'{self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()}/ directory was created.')
                subprocess.run(f'cd {self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()}', shell=True, text=True)
                self.DoGenInputFile()
                self.DoOptionInputFile()
                print('Inputs files was created.')
                self.TabStart.BtnCreate.setEnabled(False)
                subprocess.run(f'', shell=True, text=True)

                # if self.TabStart.CheckOrder.CheckParam.isChecked():
                #     command = 'cd '+self.TabSimuSets.SimuPath.EditPath.text()+self.TabSimuSets.SimuName.EditParam.text()+';chmod u+x '+self.TabSimuSets.InputFileName.EditParam.text()+';'+self.TabStart.StartOrder.EditParam.text()
                #     print(command)
                #     result = subprocess.run(command, shell=True, text=True)
                #     print('Simulation launched')
                #     error = result.stderr
                #     if type(error)!= type(None):
                #         print(result.stderr)
                #         print('Simulation not launched but you can still launch yourself the input shell file created in the desired directory.\n')
                # else:
                #     print('All you have to do is launch the input shell file created in the desired directory.')

                    
    


    def DoGenInputFile(self):
        with open(self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()+'/'+self.TabSimuFiles.InGenFileName.EditParam.text(), "w") as file:
            file.write(f'{self.DirPath}/../../Algorithm/gen_tout_multi <<!') # Header
            file.write('\n')
            file.write(self.TabSimuSets.Algo.ComboParam.currentText())
            if self.TabStart.CheckParallel.CheckParam.isChecked(): file.write('_parp')
            file.write(' # Algorithm')
            file.write('\n')
            if self.TabStart.CheckParallel.CheckParam.isChecked(): file.write(str(self.NbCoresPerSimu))
            else: file.write(str(self.TabStart.NbCores.SpinParam.value()))
            file.write(' # Number of cores')
            file.write('\n')
            file.write(str(self.TabOrbitsParams.Units.ComboParam.currentIndex()))
            file.write(f' # Units: {self.TabOrbitsParams.Units.ComboParam.currentText()}')
            file.write('\n')
            file.write(str(self.TabOrbitsParams.Coord.ComboParam.currentIndex()))
            file.write(f' # Coordinate: {self.TabOrbitsParams.Coord.ComboParam.currentText()}')
            file.write('\n')
            if self.TabSimuSets.CheckPartRemovedBody.CheckParam.isChecked(): 
                file.write('0')
                file.write(' # No radius for bodies other than the central body')
            else: 
                file.write(str(self.TabSimuSets.RminBodyType.ComboParam.currentIndex()+1))
                file.write(f' # {self.TabSimuSets.RminBodyType.ComboParam.currentText()} for bodies other than the central body')
            file.write('\n')
            if self.TabSimuSets.RminBodyType.ComboParam.currentIndex()==1:
                file.write(str(self.TabSimuSets.RminBody.SpinParam.value()))
                file.write(' # Multiple of the Hill radius')
                file.write('\n')
            file.write(str(round(float(self.TabOrbitsParams.TablePriors.item(0,0).text())*0.000954588, 3)))
            file.write(' # Mass of central body [Msun]')
            file.write('\n')
            file.write(str(self.TabOrbitsParams.NbBodies.SpinParam.value()-1))
            file.write(' # Number of bodies')
            file.write('\n')
            for i in range(1, self.TabOrbitsParams.NbBodies.SpinParam.value()):
                for j in range(len(self.TabOrbitsParams.LabelParams)):
                    file.write(self.TabOrbitsParams.TablePriors.item(i, j).text())
                    file.write(' ')
                file.write(f' # Initial orbit {i} parameters (m[Mjup] a[AU] e i[°] Om[°] om[°] M)')
                file.write('\n')
            file.write(self.TabSimuFiles.InBodFileName.EditParam.text())
            file.write(' # Input bodies file')
            file.write('\n')
            file.write(str(self.TabOrbitsParams.RandSeed.SpinParam.value()))
            file.write(' # Random seed')
            file.write('\n')
            file.write(self.TabSimuFiles.InPartFileName.EditParam.text())
            file.write(' # Input particules file')
            file.write('\n')
            file.write(self.TabSimuFiles.InSetFileName.EditParam.text())
            file.write(' # Input settings file')
            file.write('\n')
            file.write(str(self.TabOrbitsParams.NbPart.SpinParam.value()))
            file.write(' # Number of particles')
            file.write('\n')
            if self.TabOrbitsParams.NbPart.SpinParam.value()!=0:
                if self.TabStart.CheckParallel.CheckParam.isChecked(): file.write(str(self.NbPartPerSimu))
                else: file.write(str(self.TabOrbitsParams.NbPart.SpinParam.value()))
                file.write(' # Number of particles per simulation')
                file.write('\n')
                file.write(self.TabOrbitsParams.eMin.SpinParam.text()+' '+self.TabOrbitsParams.eMax.SpinParam.text())
                file.write(' # Range of particules eccentricity')
                file.write('\n')
                file.write(self.TabOrbitsParams.iMax.SpinParam.text())
                file.write(' # Maximum of inclinaison [°]')
                file.write('\n')
                file.write(self.TabOrbitsParams.aMin.SpinParam.text()+' '+self.TabOrbitsParams.aMax.SpinParam.text())
                file.write(' # Range of particles half major axis [AU]')
                file.write('\n')
                file.write('!')



    def DoOptionInputFile(self):
        with open(self.TabSimuFiles.SimuPath.EditPath.text()+self.TabSimuFiles.SimuName.EditParam.text()+'/'+self.TabSimuFiles.InSetFileName.EditParam.text(), "w") as file:
            file.write(str(self.TabSimuSets.T0.SpinParam.value()))
            file.write(' ')
            file.write(str(self.TabSimuSets.TMax.SpinParam.value()))
            file.write(' ')
            file.write(str(self.TabSimuSets.dT.SpinParam.value()))
            file.write(' # Initial, final and step time [yr]')
            file.write('\n')
            file.write(str(self.TabSimuFiles.OutFreq.SpinParam.value()))
            file.write(' ')
            file.write(str(self.TabSimuFiles.DumpFreq.SpinParam.value()))
            file.write(' # Output and dump frquencies [yr]')
            file.write('\n')
            if self.TabSimuSets.CheckSpheBody0.CheckParam.isChecked(): file.write('T')
            else: file.write('F')
            file.write(' ')
            if self.TabSimuSets.CheckPartRemovedBody.CheckParam.isChecked() or self.TabSimuSets.CheckPartRemovedBody0.CheckParam.isChecked(): file.write('T')
            else: file.write('F')
            file.write(' ')
            if self.TabSimuSets.CheckJacobiInt.CheckParam.isChecked(): file.write('T')
            else: file.write('F')
            file.write(' ')
            if self.TabSimuSets.CheckEandL.CheckParam.isChecked(): file.write('T')
            else: file.write('F')
            file.write(' ')
            if self.TabSimuSets.OutputDataType.ComboParam.currentIndex()==0: file.write('T F')
            else : file.write('F T')
            file.write(' # Options')
            file.write('\n')
            if self.TabSimuSets.CheckPartRemovedBody.CheckParam.isChecked() or self.TabSimuSets.CheckPartRemovedBody0.CheckParam.isChecked():
                if self.TabSimuSets.CheckPartRemovedBody0.CheckParam.isChecked(): 
                    file.write(str(int(self.TabSimuSets.RminBody0.SpinParam.value())))
                    file.write(' ')
                    file.write(str(int(self.TabSimuSets.RmaxBody0.SpinParam.value())))
                    file.write(' ')
                    file.write(str(int(self.TabSimuSets.RmaxBody0Unbound.SpinParam.value())))
                    file.write(' ')
                    file.write(str(int(self.TabSimuSets.PminBody0.SpinParam.value())))
                else : file.write('-1 -1 -1 -1')
                if self.TabSimuSets.CheckPartRemovedBody.CheckParam.isChecked(): file.write('T')
                else: file.write('F')
                file.write(' # Parameters for removing particles that are too far from or too close to bodies [AU]')
                file.write('\n')
            file.write(self.TabSimuFiles.SimuPath.EditPath.text())
            file.write(' # Parent path of the simulation directory')
            file.write('\n')
            file.write(self.TabSimuFiles.SimuName.EditParam.text())
            file.write(' # Name of the simulation directory')
            file.write('\n')
            file.write(self.TabSimuFiles.OutputFileName.EditParam.text())
            file.write('\n')
            file.write('new')





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetNewSimu()
    WindowParam.show()
    app.exec() # Application execution
