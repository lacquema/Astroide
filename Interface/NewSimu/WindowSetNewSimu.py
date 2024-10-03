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
        self.Tab1 = TabSimuFiles()
        self.Container.addTab(self.Tab1, 'Simulation files')
        self.Tab1.SimuPath.EditPath.textChanged.connect(self.ChangeStartOrder)
        self.Tab1.SimuName.EditParam.textChanged.connect(self.ChangeStartOrder)

        # Tab 2
        self.Tab2 = TabSimuSet()
        self.Container.addTab(self.Tab2, 'Simulation settings')

        # Tab 3
        self.Tab3 = TabOrbitParamSet()
        self.Container.addTab(self.Tab3, 'Orbital parameters')
        self.Tab3.NbPart.SpinParam.valueChanged.connect(self.ChangeSumaPara)

        # Tab 4
        self.Tab4 = TabStartSet()
        self.Tab4.BtnStart.clicked.connect(self.StartSimulation)
        self.Container.addTab(self.Tab4, 'Starting')
        self.Tab4.NbHours.SpinParam.valueChanged.connect(self.ChangeStartOrder)
        self.Tab4.NbCores.SpinParam.valueChanged.connect(self.ChangeStartOrder)
        self.Tab4.NbCores.SpinParam.valueChanged.connect(self.ChangeSumaPara)
        self.Tab4.NbBranch.SpinParam.valueChanged.connect(self.ChangeSumaPara)
        self.ChangeSumaPara()
        self.Tab4.BtnReset.clicked.connect(self.ChangeSumaPara)
        

        # Container
        self.setCentralWidget(self.Container)

        # Status bar
        self.setStatusBar(QStatusBar(self))


    def ChangeSumaPara(self):
        self.NbCoresPerSimu = self.Tab4.NbCores.SpinParam.value()//self.Tab4.NbBranch.SpinParam.value()
        self.NbPartPerSimu = int(np.ceil(self.Tab3.NbPart.SpinParam.value()/self.Tab4.NbBranch.SpinParam.value()))
        self.Tab4.SumaPara.setText(f'=> {self.NbPartPerSimu} particules and {self.NbCoresPerSimu} cores per simulation')


    # Emition of the CloseEvent signal when the parameter window is closed
    def closeEvent(self, e):
        self.SignalCloseWindowSetNewSimu.emit() 

    def ChangeStartOrder(self):
        self.Tab4.StartOrder.EditParam.setText(f'oarsub -l nodes=1/core={self.NbCoresPerSimu},walltime={self.Tab4.NbHours.SpinParam.value()} --project dynapla {self.Tab1.SimuPath.EditPath.text()+self.Tab1.SimuName.EditParam.text()}/{self.Tab1.StartFileName.EditParam.text()}')
    
    def StartSimulation(self):
        print('--------------------------')
        if len(self.Tab2.SimuPath.EditPath.text())==0:
            print('Simulation path not given.')
            print('Check your inputs.')
        else:
            if os.path.exists(self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()):
                print('This directory already exists.')
                print('Check your inputs.')
            else: 
                os.makedirs(self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text())
                if len(self.Tab3.PathData.EditPath.text())==0:
                    os.rmdir(self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text())
                    print('Data file not given.')
                    print('Check your inputs.')
                else:
                    if self.Tab3.FormatAstro.ComboParam.currentIndex() == 0:
                        os.rmdir(self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text())
                        print('Astrometric data format not given.')
                        print('Check your inputs.')
                    else:
                        print(f'{self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()}/ directory was created.')
                        subprocess.run(f'cd {self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()}', shell=True, text=True)
                        shutil.copy(self.Tab3.PathData.EditPath.text(), self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()+'/'+self.Tab3.DataFileName.EditParam.text())                    
                        print('Data file was copied.')
                        self.DoInputShell()
                        print('Input shell file was created.')
                        if self.Tab4.CheckOrder.CheckParam.isChecked():
                            command = 'cd '+self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()+';chmod u+x '+self.Tab2.InputFileName.EditParam.text()+';'+self.Tab4.StartOrder.EditParam.text()
                            print(command)
                            result = subprocess.run(command, shell=True, text=True)
                            print('Simulation launched')
                            error = result.stderr
                            if type(error)!= type(None):
                                print(result.stderr)
                                print('Simulation not launched but you can still launch yourself the input shell file created in the desired directory.\n')
                        else:
                            print('All you have to do is launch the input shell file created in the desired directory.')

                    
    


    def DoInputShell(self):
        with open(self.Tab2.SimuPath.EditPath.text()+self.Tab2.SimuName.EditParam.text()+'/'+self.Tab2.InputFileName.EditParam.text(), "w") as file:
            file.write('#! /bin/bash\nexport OMP_NUM_THREADS=8\nexport STACKSIZE=1000000\n'+self.DirPath+'/../../Algorithm/bin/astrom_mcmcop <<!') # Header
            file.write('\n')
            if self.Tab4.CheckLM.isChecked() and self.Tab4.CheckMCMC.isChecked(): # Choice
                file.write('2')
                file.write(' # New simulation LM and MCMC')
            elif not self.Tab4.CheckLM.isChecked() and self.Tab4.CheckMCMC.isChecked():
                file.write('3')
                file.write(' # New simulation LM only')
            elif self.Tab4.CheckLM.isChecked() and not self.Tab4.CheckMCMC.isChecked():
                file.write('4')
                file.write(' # New simulation MCMC only')
            file.write('\n')
            if self.Tab3.RelAstro.CheckParam.isChecked(): file.write('1 ') # Data
            else: file.write('0 ')
            if self.Tab3.AbsRV.CheckParam.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.Tab3.RelRV.CheckParam.isChecked(): file.write('1 ')
            else: file.write('0 ')
            if self.Tab3.AbsAstro.CheckParam.isChecked(): file.write('1')
            else: file.write('0')
            file.write(' # Type of data (RelAstro? AbsRV? RelRV? AbsAstro?)')
            file.write('\n')
            file.write(str(self.Tab3.NbOrbitsValue)) # Number of orbits
            file.write(' # Number of orbits')
            file.write('\n')
            if self.Tab3.AbsAstro.CheckParam.isChecked() or self.Tab3.RelAstro.CheckParam.isChecked() or self.Tab3.AbsRV.CheckParam.isChecked() or self.Tab3.RelRV.CheckParam.isChecked(): 
                file.write(self.Tab3.DataFileName.EditParam.text())
                # file.write(' # Data file')
                file.write('\n')
            file.write('1d-'+str(self.Tab2.Precision.SpinParam.value())) # Precision of simulation
            file.write(' # Precision')
            file.write('\n')
            file.write(str(self.Tab3.FormatDate.ComboParam.currentIndex()+1))
            file.write(' '+str(self.Tab3.FormatAstro.ComboParam.currentIndex()))
            if self.Tab3.CheckCorrCoef.CheckParam.isChecked(): file.write(' 1')
            else: file.write(' 0')
            file.write(' # Format data (1:DDMMYYYY/2:JD 1:(DEC,RA)/2:(SEP,PA) CorrCoeff?')
            file.write('\n')
            if self.Tab3.RelRV.CheckParam.isChecked() or self.Tab3.AbsRV.CheckParam.isChecked():
                if self.Tab3.CheckJitter.isChecked(): 
                    file.write('1')
                else: 
                    file.write('0')
                file.write(' # Jitter?')
                file.write('\n')
            file.write(self.Tab3.SystDist.SpinParam.text()+' '+self.Tab3.SystDistUnit.ComboParam.currentText())
            file.write(' # Distance')
            file.write('\n')
            file.write(str(float(self.Tab3.TablePriors.item(0,0).text())*0.000954588))
            file.write(' # First guess of center mass (ms)')
            file.write('\n')
            file.write(self.Tab2.OutFileName.EditParam.text())
            # file.write(' # Result file')
            file.write('\n')
            file.write(self.Tab2.DumpFileName.EditParam.text())
            # file.write(' # Dump file')
            file.write('\n')
            file.write(self.Tab2.DumpFreq.SpinParam.text())
            file.write(' # Dump frequency')
            file.write('\n')
            c=0 # Number of mass prior
            for i in range(len(self.Tab3.ListPriorMass)):
                DistribIndex = self.Tab3.ListPriorMass[i].Layout.itemAt(3*self.Tab3.NbBodies.SpinParam.value()+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0: c+=1
            file.write(str(c))
            file.write(' # Number of masses prior')
            file.write('\n')
            for i in range(len(self.Tab3.ListPriorMass)):
                DistribIndex = self.Tab3.ListPriorMass[i].Layout.itemAt(3*self.Tab3.NbBodies.SpinParam.value()+1).widget().ComboParam.currentIndex()
                if DistribIndex != 0:
                    for j in range(self.Tab3.NbBodies.SpinParam.value()):
                        file.write(self.Tab3.ListPriorMass[i].Layout.itemAt(3*j+1).widget().SpinParam.text()+' ')
                    file.write(' # Coefficients')
                    file.write('\n')
                    file.write(str(DistribIndex))
                    file.write(' # Distribution (1:Normal, 2:Log, 3:Uniform, 4:Fixed)')
                    file.write('\n')
                    if DistribIndex == 1 or DistribIndex == 2:
                        file.write(self.Tab3.ListPriorMass[i].Mean.SpinParam.text())
                        file.write(' ')
                        file.write(self.Tab3.ListPriorMass[i].SD.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 3:
                        file.write(self.Tab3.ListPriorMass[i].Min.SpinParam.text())
                        file.write(' ')
                        file.write(self.Tab3.ListPriorMass[i].Max.SpinParam.text())
                        file.write(' ')
                    if DistribIndex == 4:
                        file.write(self.Tab3.ListPriorMass[i].Value.SpinParam.text())
                        file.write(' ')
                    file.write(self.Tab3.ListPriorMass[i].PriorUnit.ComboParam.currentText())
                    file.write(' # Distribution parameters')
                    file.write('\n')
            file.write(self.Tab3.RefTime.SpinParam.text())
            file.write(' # Reference of time')
            file.write('\n')
            if self.Tab3.CheckJitter.isChecked(): 
                file.write(self.Tab3.Jitter.SpinParam.text()+' '+self.Tab3.V0.SpinParam.text())
                file.write(' # Initial VO and Jitter')
                file.write('\n')
            for i in range(1, self.Tab3.NbBodies.SpinParam.value()):
                for j in range(len(self.Tab3.LabelParams)):
                    file.write(self.Tab3.TablePriors.item(i, j).text())
                    if j == 0: file.write(' mj')
                    file.write(' ')
                file.write(' # First guess of orbit parameters (m[mj] a[AU] e i[deg] Om[deg] om[deg] tp[MJD])')
                file.write('\n')
            file.write(self.Tab3.PMin.SpinParam.text()+' '+self.Tab3.PMax.SpinParam.text())
            file.write(' # Range of permited period')
            file.write('\n')
            file.write(self.Tab3.aMin.SpinParam.text()+' '+self.Tab3.aMax.SpinParam.text())
            file.write(' # Range of permited half major axis')
            file.write('\n')
            file.write(self.Tab3.eMin.SpinParam.text()+' '+self.Tab3.eMax.SpinParam.text())
            file.write(' # Range of permited eccentricity')
            file.write('\n')
            file.write('exit')
            file.write('\n')
            file.write('!')
            
            



            


            





# Check
if __name__=="__main__":
    app = QApplication(sys.argv) # Application creation
    WindowParam = WindowSetNewSimu()
    WindowParam.show()
    app.exec() # Application execution
