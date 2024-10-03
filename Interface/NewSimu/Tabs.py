#! /Users/lacquema/ByeGildas/bin/python3


### --- Packages --- ###

# Transverse packages
import sys
import os

# PyQt packages
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtWidgets import QTableWidgetItem, QTableWidget, QFileDialog, QTabWidget, QScrollArea, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QWidget, QStatusBar, QApplication, QProgressBar, QLabel, QCheckBox
from PyQt6.QtGui import QIcon, QFont

# My packages
from Parameters import *
from UtilsNewSimu import DelAllWidgetsBtw




class GeneralTab(QWidget):
    
    def __init__(self):
        super().__init__()

        # Layout initialisation
        self.Layout = QVBoxLayout()

        # Reset button
        self.BtnReset = QPushButton('Reset')
        self.BtnReset.setStatusTip('Reset tab settings')
        self.Layout.addWidget(self.BtnReset)
        self.BtnReset.clicked.connect(self.ResetParams)
        
        self.InitWidgets()

        self.Layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        
        self.setLayout(self.Layout)

    # Widgets initialisations
    def InitWidgets(self):
        return
    
    # Reset all widgets of the parameters window
    def ResetParams(self):
        DelAllWidgetsBtw(self.Layout, 1, self.Layout.count())
        self.InitWidgets()


    def ValidatedIfIn(self, WidgetEditing, ListValidCharacters):
        TextNew = WidgetEditing.text()
        c = 0
        for i in range(len(TextNew)):
            if TextNew[i] not in  ListValidCharacters:
                c=+1
                break
        if c != 0:
            WidgetEditing.setText(self.TextOld)
            WidgetEditing.setCursorPosition(i)
        else:
            self.TextOld = TextNew
    


class TabSimuFiles(GeneralTab):
    
    def __init__(self):
        super().__init__()

    def InitWidgets(self):

        
        self.SimuPath = PathBrowser('Path', 'Path where create the adjustment directory', 0)
        self.Layout.addWidget(self.SimuPath)

        self.SimuPath.Layout.addSpacing(20)

        self.SimuName = LineEdit('Directory', 'Name you want to give to the adjustment directory', '')
        self.SimuPath.Layout.addWidget(self.SimuName)

        self.Layout.addWidget(Delimiter(Title='Input files:'))

        self.StartFileName = LineEdit('General input file', 'Name you want to give to the general entry which is the shell startup file with the extension', 'start.sh')
        self.Layout.addWidget(self.StartFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.InBodFileName = LineEdit('Bodies input file', 'Name you want to give to the input file for bodies with the extension', 'bodies.in')
        self.Layout.addWidget(self.InBodFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.InParFileName = LineEdit('Particules input file', 'Name you want to give to the input file for test particules with the extension', 'particules.in')
        self.Layout.addWidget(self.InParFileName, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Delimiter(Title='Output files:'))

        self.OutputFileName = LineEdit('Output binary file', 'Name you want to give to the output binary file (no extension)', 'output')
        self.Layout.addWidget(self.OutputFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        self.OutputFileName.Layout.addSpacing(100)
        self.DumpFreq = SpinBox('Frequency', 'Time between two saves to the output file [yr]', 10000, 1, None, 1)
        self.OutputFileName.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        self.DumpFileName = LineEdit('Dump binary file', 'Name you want to give to the dump binary file (no extension)', 'dump')
        self.Layout.addWidget(self.DumpFileName, alignment=Qt.AlignmentFlag.AlignLeft)
        self.DumpFileName.Layout.addSpacing(100)
        self.DumpFreq = SpinBox('Frequency', 'Time between two saves to the dump file [yr]', 100000, 1, None, 1)
        self.DumpFileName.Layout.addWidget(self.DumpFreq, alignment=Qt.AlignmentFlag.AlignLeft)

        self.Layout.setSpacing(0)


class TabSimuSet(GeneralTab):
    
    def __init__(self):
        super().__init__()

    def InitWidgets(self):

        self.SimuType = ComboBox('Algorithm', 'Algorithm to be used for integration', ['whm', 'whm_s6b', 'rmvs3'])
        self.SimuType.ComboParam.setCurrentIndex(2)
        self.Layout.addWidget(self.SimuType)

        self.T0 = DoubleSpinBox('Integration times', 'Time when integration begins [yr]', 0, 0)
        self.Layout.addWidget(self.T0)
        self.T0.Layout.addWidget(QLabel('     <-'))
        self.dT = DoubleSpinBox(None, 'Time between two integration steps [yr] (0.05 for 20% of the shortest orbit)', 0.05)
        self.T0.Layout.addWidget(self.dT)
        self.T0.Layout.addWidget(QLabel('->     '))
        self.TMax = DoubleSpinBox(None, 'Time when integration stops [yr]', 20000000, 0)
        self.T0.Layout.addWidget(self.TMax)

        self.OutputDataType = ComboBox('Output data type', 'Type in which output data is written', ['REAL*4', 'INTEGER*2'])
        self.Layout.addWidget(self.OutputDataType)

        self.Layout.addWidget(Delimiter(Title='Options:'))

        self.CheckSpheStar = CheckBox('Central body sphericity', 'Sphericity of the central body')
        self.CheckSpheStar.CheckParam.stateChanged.connect(self.ChangeStateSpheStar)
        self.Layout.addWidget(self.CheckSpheStar)
        self.CheckSpheStar.Layout.addSpacing(60)
        self.J2 = DoubleSpinBox('J2', 'Quadrupole graviational moment')
        self.J2.setEnabled(self.CheckSpheStar.CheckParam.isChecked())
        self.CheckSpheStar.Layout.addWidget(self.J2)
        self.CheckSpheStar.Layout.addSpacing(20)
        self.J4 = DoubleSpinBox('J4', 'Octopole graviational moment')
        self.J4.setEnabled(self.CheckSpheStar.CheckParam.isChecked())
        self.CheckSpheStar.Layout.addWidget(self.J4)
        self.CheckSpheStar.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        self.Layout.addWidget(Label('Particle removal:'))
        self.CheckPartRemovedBody0 = CheckBox('by the central body', 'Remove particles that are too close or too far to the central body')
        self.CheckPartRemovedBody0.CheckParam.stateChanged.connect(self.ChangeStatePartRemovedBody0)
        self.Layout.addWidget(self.CheckPartRemovedBody0)
        self.CheckPartRemovedBody0.Layout.insertWidget(0, Label('   '), alignment=Qt.AlignmentFlag.AlignLeft)
        self.CheckPartRemovedBody0.Layout.addSpacing(60)
        self.RminBody0 = DoubleSpinBox('Rmin', 'Distance below which a particle is removed because it is too close to the central body (-1 to ignore)', -1)
        self.RminBody0.setEnabled(self.CheckPartRemovedBody0.CheckParam.isChecked())
        self.CheckPartRemovedBody0.Layout.addWidget(self.RminBody0)
        self.CheckPartRemovedBody0.Layout.addSpacing(20)
        self.RmaxBody0 = DoubleSpinBox('Rmax', 'Distance beyond which a particle with a bound orbit is removed because it is too far to the central body (-1 to ignore)', -1)
        self.RmaxBody0.setEnabled(self.CheckPartRemovedBody0.CheckParam.isChecked())
        self.CheckPartRemovedBody0.Layout.addWidget(self.RmaxBody0)
        self.RmaxBody0Unbound = DoubleSpinBox(None, 'Distance beyond which a particle with an unbound orbit is removed because it is too far from the central body (-1 to ignore)', -1)
        self.RmaxBody0Unbound.setEnabled(self.CheckPartRemovedBody0.CheckParam.isChecked())
        self.CheckPartRemovedBody0.Layout.addWidget(self.RmaxBody0Unbound)
        self.CheckPartRemovedBody0.Layout.addSpacing(20)
        self.PminBody0 = DoubleSpinBox('Pmin', 'Periastron below which a particle is removed because it is too close to the central body (-1 to ignore)', -1)
        self.PminBody0.setEnabled(self.CheckPartRemovedBody0.CheckParam.isChecked())
        self.CheckPartRemovedBody0.Layout.addWidget(self.PminBody0)
        self.CheckPartRemovedBody0.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        self.CheckPartRemovedBody = CheckBox('by the other bodies', 'Remove particles that are too close to a body')
        self.CheckPartRemovedBody.CheckParam.stateChanged.connect(self.ChangeStatePartRemovedBody)
        self.Layout.addWidget(self.CheckPartRemovedBody)
        self.CheckPartRemovedBody.Layout.insertWidget(0, Label('   '), alignment=Qt.AlignmentFlag.AlignLeft)
        self.CheckPartRemovedBody.Layout.addSpacing(60)
        self.DistType = ComboBox(None, 'Type of distance', ['Physical radius', 'Multiple of Hill radius'])
        self.DistType.setEnabled(self.CheckPartRemovedBody.CheckParam.isChecked())
        self.CheckPartRemovedBody.Layout.addWidget(self.DistType)
        self.CheckPartRemovedBody.Layout.addSpacing(20)
        self.RminBody = DoubleSpinBox(None, 'Distance below which a particle is removed because it is too close to a body (-1 to ignore)')
        self.RminBody.setEnabled(self.CheckPartRemovedBody.CheckParam.isChecked())
        self.CheckPartRemovedBody.Layout.addWidget(self.RminBody)
        self.CheckPartRemovedBody.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        self.RBodies = ComboBox('Radius of bodies', 'Choice of the radius of bodies', ['No', 'Physical', 'Hill multiple'])

        self.AddCalculation = Label('Additional calculation:')
        self.Layout.addWidget(self.AddCalculation)
        self.AddCalculation.Layout.addSpacing(60)
        self.CheckJacobiInt = CheckBox('Jacobi integral', 'Calculate the jacobi integral and save it in a jacobi.out file')
        self.AddCalculation.Layout.addWidget(self.CheckJacobiInt)
        self.AddCalculation.Layout.addSpacing(20)
        self.CheckEandL = CheckBox('Energy and angular momentum', 'Calculate the energy and angular momentum and save them in an energy.out file')
        self.AddCalculation.Layout.addWidget(self.CheckEandL)
        self.AddCalculation.Layout.setAlignment(Qt.AlignmentFlag.AlignLeft)

        self.Layout.setSpacing(0)


    def ChangeStateSpheStar(self, state):
        self.J2.setEnabled(state)
        self.J4.setEnabled(state)

    def ChangeStatePartRemovedBody0(self, state):
        self.RminBody0.setEnabled(state)
        self.RmaxBody0.setEnabled(state)
        self.RmaxBody0Unbound.setEnabled(state)
        self.PminBody0.setEnabled(state)

    def ChangeStatePartRemovedBody(self, state):
        self.DistType.setEnabled(state)
        self.RminBody.setEnabled(state)





class TabOrbitParamSet(GeneralTab):
    
    def __init__(self):
        super().__init__()


    def InitWidgets(self):

        self.WidgetH = QWidget()
        self.LayoutH = QHBoxLayout()
        
        self.LayoutV1 = QVBoxLayout()

        self.NbBodies = SpinBox('Number of bodies', 'Number of bodies', 2, 2, 5, 1)
        self.LayoutV1.addWidget(self.NbBodies, alignment=Qt.AlignmentFlag.AlignLeft)
        self.NbBodiesValue = self.NbBodies.SpinParam.value()
        self.NbBodies.SpinParam.valueChanged.connect(self.ChangeNbBodies)

        self.NbOrbitsValue = self.NbBodiesValue-1
        self.NbOrbits = QLabel(f'=> {self.NbOrbitsValue} Orbit')
        self.NbBodies.Layout.addWidget(self.NbOrbits)

        self.Units = ComboBox('System units', 'Choice of the units used', ["AU, Msun, yr/2π", 'AU, Msun/4π^2, yr'])
        self.LayoutV1.addWidget(self.Units)

        self.Coord = ComboBox('Coordinate', 'Choice of the coordinate system used', ['Ecliptic', 'Invariable plane'])
        self.LayoutV1.addWidget(self.Coord)

        self.LayoutV1.addWidget(Delimiter(Title='Initial orbit parameters of particles :'))

        self.NbPart = SpinBox('Number of particles', 'Number of particles at start', 100000, 0, None)
        self.LayoutV1.addWidget(self.NbPart)



        self.RandSeed = SpinBox('Random seed', 'Seed used to create a random value', 654876543, 1, None)
        self.LayoutV1.addWidget(self.RandSeed)

        self.aMin = DoubleSpinBox('Semi-major axis', 'Mimimum of orbits semi-major axis [AU]', 15, 0, None)
        self.LayoutV1.addWidget(self.aMin)
        self.aMin.Layout.addWidget(QLabel('   <->'))
        self.aMax = DoubleSpinBox(None, 'Maximum of orbits semi-major axis [AU]', 80, 0, None)
        self.aMin.Layout.addWidget(self.aMax)

        self.eMin = DoubleSpinBox('Eccentricity', 'Mimimum of orbits eccentricity', 0, 0, 1)
        self.LayoutV1.addWidget(self.eMin)
        self.eMin.Layout.addWidget(QLabel('   <->'))
        self.eMax = DoubleSpinBox(None, 'Maximum of orbits eccentricity', 0.05, 0, 1)
        self.eMin.Layout.addWidget(self.eMax)

        self.iMax = DoubleSpinBox('Inclinaison max', 'Maximum of orbits inclinaison [°]', 2, 0, 360)
        self.LayoutV1.addWidget(self.iMax)

        self.LayoutV1.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.LayoutV1.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.LayoutV1.setSpacing(0)
        self.LayoutH.addLayout(self.LayoutV1)

        self.LayoutH.addSpacing(20)

        self.LayoutV2 = QVBoxLayout()

        self.LayoutV2.addWidget(Delimiter(Title = 'Initial orbit parameters of each bodies:', Width=0))

        self.TablePriors = QTableWidget()
        self.TablePriors.setStatusTip('First guess of orbits parameters of each bodies.')
        self.TablePriors.setRowCount(self.NbBodiesValue)
        self.LabelParams = ['m [Mjup]', 'a [AU]', 'e', 'i [°]', 'w [°]', 'W [°]', 'M']
        self.TablePriors.setColumnCount(len(self.LabelParams))
        self.TablePriors.setHorizontalHeaderLabels(self.LabelParams)
        self.LayoutV2.addWidget(self.TablePriors, alignment=Qt.AlignmentFlag.AlignVCenter)
        for i in range(self.NbBodiesValue):
            for j in range(len(self.LabelParams)):
                if i==0 and j!=0: self.TablePriors.setItem(i, j, QTableWidgetItem('X'))
                else: self.TablePriors.setItem(i, j, QTableWidgetItem('0.'))
                self.TablePriors.item(i, j).setTextAlignment(Qt.AlignmentFlag.AlignCenter)

        self.TablePriors.itemChanged.connect(self.ValidationItemTbl)
        self.TablePriors.cellClicked.connect(self.SaveOldTextTbl)

        # self.RefTime = SpinBox('Time reference', 'Reference of the time to count the orbit phase [MJD]', 0, 0, None, 1)
        # self.LayoutV2.addWidget(self.RefTime, alignment=Qt.AlignmentFlag.AlignLeft)

        self.LayoutV2.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.LayoutH.addLayout(self.LayoutV2)

        self.WidgetH.setLayout(self.LayoutH)
        self.Layout.addWidget(self.WidgetH)

        self.Layout.setSpacing(0)



    def ChangeNbBodies(self):
        # Change widgets
        OldValue = self.NbBodiesValue
        self.NbBodiesValue = self.NbBodies.SpinParam.value()
        self.NbOrbitsValue = self.NbBodiesValue-1
        if  self.NbOrbitsValue==1: self.NbOrbits.setText(f'=> {self.NbOrbitsValue} Orbit')
        else: self.NbOrbits.setText(f'=> {self.NbOrbitsValue} Orbits')

        # Change the number of row in the priors table
        self.TablePriors.setRowCount(self.NbBodiesValue)

        if self.NbBodiesValue > OldValue:
            for n in range(len(self.LabelParams)):
                self.TablePriors.setItem(self.NbBodiesValue-1, n, QTableWidgetItem('0.'))
                self.TablePriors.item(self.NbBodiesValue-1, n).setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            
            
    def ValidationItemTbl(self):
        if self.TablePriors.currentRow() != -1 and self.TablePriors.currentColumn() != -1: # corresponding of initial state (-1,-1) of the table where item is None
            TextNew = self.TablePriors.currentItem().text()
            ListValidCharacters = '1234567890.'
            PointSuppr = False
            if len(TextNew)==0:
                self.TablePriors.currentItem().setText(self.TextOld)
            else:
                for i in range(len(TextNew)):
                    if TextNew[i] not in  ListValidCharacters:
                        self.TablePriors.currentItem().setText(self.TextOld)
                        break
                    if TextNew[i]=='.' and PointSuppr==False: 
                            ListValidCharacters='1234567890'
                            PointSuppr = True 
                if self.TablePriors.currentRow() == 0 and self.TablePriors.currentColumn() != 0:
                    self.TablePriors.currentItem().setText(self.TextOld)

    def SaveOldTextTbl(self):
        self.TextOld = self.TablePriors.currentItem().text()



class TabStartSet(GeneralTab):
    
    def __init__(self):
        super().__init__()


    def InitWidgets(self):
    
        self.NbCores = SpinBox('Number of cores', 'Number of cores to be used', 48, 1, None, 1)
        self.Layout.addWidget(self.NbCores, alignment=Qt.AlignmentFlag.AlignLeft)

        self.CheckParallel = CheckBox('Parallelization', 'Parallelization of the simulation algorithm')
        self.Layout.addWidget(self.CheckParallel)
        self.CheckParallel.Layout.addSpacing(60)
        self.NbBranch = SpinBox('Number of branches', 'Number of parallelization branches', 6, 1, None)
        self.CheckParallel.Layout.addWidget(self.NbBranch)
        self.CheckParallel.CheckParam.stateChanged.connect(self.CheckParallelChange)
        self.CheckParallel.Layout.addSpacing(20)
        self.SumaPara = QLabel('')
        self.CheckParallel.Layout.addWidget(self.SumaPara)


        self.BtnCreate = QPushButton('Create startup files')
        self.Layout.addWidget(self.BtnCreate, alignment=Qt.AlignmentFlag.AlignRight)

        self.CheckOrder = CheckBox('Starting order', 'If you just want to create the input file, but dont want to run the command in the terminal')
        self.Layout.addWidget(self.CheckOrder)
        self.CheckOrder.CheckParam.stateChanged.connect(self.CheckStartOrderChange)
        
        self.NbHours = SpinBox('Simulation duration', 'Simulation duration [hour]', 48, 1, None)
        self.Layout.addWidget(self.NbHours, alignment=Qt.AlignmentFlag.AlignLeft)

        self.StartOrderValue = f'oarsub -l nodes=1/core=8,walltime={self.NbHours.SpinParam.value()} --project dynapla ./input.sh'
        self.StartOrder = LineEdit('Order', 'Terminal order to start the adjustment', self.StartOrderValue)
        self.Layout.addWidget(self.StartOrder)
        
        self.BtnStart = QPushButton('Start the simulation')
        self.Layout.addWidget(self.BtnStart, alignment=Qt.AlignmentFlag.AlignRight)

        self.CheckParallelChange(self.CheckParallel.CheckParam.isChecked())
        self.CheckStartOrderChange(self.CheckOrder.CheckParam.isChecked())
        

    def CheckParallelChange(self, state):
        self.NbBranch.setEnabled(state)
        self.SumaPara.setEnabled(state)


    def CheckStartOrderChange(self, state):
        self.StartOrder.setEnabled(state)
        self.BtnStart.setEnabled(state)
        self.NbHours.setEnabled(state)


        
    
