# Example script to generate sub-simulations for the HJS integrator.
# This script sets up a system with one star, two planets, and two disk regions.
# Adjust parameters below as needed for your specific simulation setup.

# Path to the executables
BINPATH='<environment_path>/Astroide/Code/bin'

# Path to the working directory
WORKPATH='...'

# Generate the sub-simulations input files
$BINPATH/gen_multi_hjs <<!
$BINPATH
hjs_par							# Integration method
0 10000000 0.05					# Simulation times
10000 100000					# Output frequencies
F T F F T F						# Computational options
0.005 200. 10. -1. F			# Removal limits
$WORKPATH
8								# Number of CPU cores per sub-simulation
1								# Units
1								# Coordinate system
0								# Planetary radius type
3								# Number of bodies
0.47							# Mass of the body 1
0,000190918						# Mass of the body 2
0,000286377						# Mass of the body 3
-1 1 0							# Orbital hierarchy of the body 2
19. 0. 0. 0. 0. 0.	    		# Initial orbit parameters of the body 2
-1 -1 1							# Orbital hierarchy of the body 3
52 0. 0. 0. 0. 0.				# Initial orbit parameters of the body 3
654876543						# Random seed
50000							# Number of debris in the disk 1
25000							# Number of debris per sub-simulation
-1 0 0							# Orbital hierarchy of the disk 1
1								# Reference frame of the disk 1
0. 0.01							# Eccentricity range in the disk in the disk 1
2.								# Maximum inclination in the disk 1
5. 52.							# Semi-major axis range in the disk 1
150000							# Number of debris in the disk 2
25000							# Number of debris per sub-simulation
-1 -1 0							# Orbital hierarchy of the disk 2
1								# Reference frame of the disk 2
0. 0.01							# Eccentricity range in the disk in the disk 2
2.								# Maximum inclination in the disk 2
52. 130.						# Semi-major axis range in the disk 2
0								# End of the inputs
!