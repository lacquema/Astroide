### This script generates the input files for the HJS integrator ###

# Path to the executables
BINPATH='/Users/lacquema/Astroide.env/Astroide/Code/bin'

# Path to the working directory
WORKPATH='/Users/lacquema/Simulations/twa7/twa7_cb_a_dyn/test'

# Generate the input files
$BINPATH/gen_multi_hjs <<!
$BINPATH
hjs								# Method (hjs / hjs_par)
0. 10000000. 0.05				# Initial, final, step time [yr]
10000. 100000.					# Output, dump frequencies [yr]
F T F F T F						# Options: int*2, real*4, nrj, debris jacobi coord, removing limits, J2 and J4 terms (T / F)
0.005 200. 10. -1. F			# Removing limits if T: rmin, rmax, rmaxu, qmin (-1=not checked), lclose (T / F)
$WORKPATH
8							    # Number of cpu cores per sub-simulation
1								# Units (0=[Msun,AU] / 1=[AU,yr])
1								# Coordinate (0=ecliptic / 1=invariant plane)
0								# Planet radius (0=None / 1=Physical / 2=Hill)
3								# Number of bodies
0.47							# Mass of the body 1 [Msun]
0.00028638						# Mass of the body 2 [Msun]
0.000324564						# Mass of the body 3 [Msun]
-1 1 0							# Orbital hierarchy of 2: body 1, body 2, body 3 (-1=inside / 0=outside / 1=body in question)
15.24 0.05 0. 0. 0. 0.			# Initial orbit parameters of 2: a[AU], e, i[deg], w[deg], W[deg], M[deg]
-1 -1 1							# Orbital hierarchy of 3: body 1, body 2, body 3 (-1=inside / 0=outside / 1=body in question)
52 0. 0. 0. 0. 0.				# Initial orbit parameters of 2: a[AU], e, i[deg], w[deg], W[deg], M[deg]
654876543						# Random seed
20000							# Number of debris in the disk 1
8000							# Number of debris per sub-simulation
-1 0 0							# Orbital hierarchy of the disk: body 1, body 2, body 3 (-1=inside / 0=outside)
1								# Coordinate of the disk (-1=orbit / 0=ecliptic / 1=centers / 2=invariant plane)
0. 0.01							# Acceptable range of eccentricity
2.								# Maximun of inclination       
5. 15.							# Acceptable range of semi-major axis
180000							# Number of debris in the disk 2
8000							# Number of debris per sub-simulation
-1 -1 0							# Orbital hierarchy of the disk: body 1, body 2, body 3 (-1=inside / 0=outside)
1								# Coordinate of the disk (-1=orbit / 0=ecliptic / 1=centers / 2=invariant plane)
0. 0.01							# Acceptable range of eccentricity
2.								# Maximun of inclination       
15. 130.						# Acceptable range of semi-major axis
0								# End of the inputs
!