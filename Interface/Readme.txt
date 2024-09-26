Credits : Version 1.1, 2024, IPAG, Antoine Lacquement
Univ. Grenoble Alpes, CNRS, IPAG, 38000 Grenoble, France

All necessary packages: sys, os, PyQt6, numpy, math, matplotlib
Watch out for the shebang beginning by '#!" if you run the program directly from the terminal.





### --- (1) Hierarchy --- ###
WindowMain <- WindowMenu, TransfertData, SnapSelector, Tools
Tools <- WindowParam, WindowPlot
WindowParam <- Parameters, Resonance, Curve





### --- (2) Routines --- ###

## --- (2.1) WindowMenu --- ##
- Necessary packages: sys, os, PyQt6
- Necessary subroutines: None
- Comment: As the data takes a while to load, a loading window will keep you waiting. 

## --- (2.2) WindowMain --- ##
- Necessary packages: sys, PyQt6
- Necessary subroutines: WindowMenu, TransfertData, SnapSelector, Tools
- Comment: This routine is the main code. This is the one to run. 
It starts by loading the data in the input zone at the top, where you must enter the different paths to the mextract.dat and followbodies.dat data files.
Meanwhile, this code displays the WindowMenu. 
It then displays the main window, in which you can change the snapshot number and plot different functions according to the parameters you've entered.

## --- (2.3) TransfertData--- ##
- Necessary packages: numpy
- Necessary subroutines: None
- Comment: It can be used to reread mextract.dat and followbodies.dat data files.

## --- (2.4) SnapSelector --- ##
- Necessary packages: sys, os, PyQt6
- Necessary subroutines: None
- Comment: This is the interactive block used to change the snapshot number.

## --- (2.5) Tools --- ##
- Necessary packages: sys, numpy, math, PyQt6
- Necessary subroutines: WindowParam, WindowPlot
- Comment: This code is the trickiest. It condenses all the functions available and links the WindowParam with the WindowPlot.
All functions have a GeneralTool parent, which provides common characteristics. And each looks for its specific parameters in the associated WindowParam.
SpaceView plots the system in 2D and 3D spaces. 
DiagramAE plots the eccentricity of the various system bodies as a function of their semi-major axis and the snapshot number.
DiagramTY plots the evolution of a planet's orbital parameter as a function of time.
DiagramXY plots the evolution of a planet's orbital parameter Y as a function of an other X.
RadProfile plots the integrated radial profile of the system particules.

## --- (2.6) WindowParam --- ##
- Necessary packages: sys, PyQt6
- Necessary subroutines: Parameters
- Comment: Window supporting the specific parameters of the selected function

## --- (2.7) WindowPlot --- ##
- Necessary packages: sys, PyQt6, matplotlib
- Necessary subroutines: None
- Comment: Window supporting the plot of the selected function

## --- (2.8) Parameters --- ##
- Necessary packages: sys, PyQt6
- Necessary subroutines: None
- Comment: All different types of parameters

## --- (2.9) Resonance --- ##
- Necessary packages: PyQt6
- Necessary subroutines: Parameters
- Comment: Initialise different resonances in the WindowParam associated

## --- (2.10) Curve --- ##
- Necessary packages: PyQt6
- Necessary subroutines: Parameters
- Comment: Initialise different curves in the WindowParam associated