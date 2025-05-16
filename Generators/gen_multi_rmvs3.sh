### This script generates the input files for the RMVS3 integrator ###

# Path to the executables
BINPATH='/Users/lacquema/Astroide.env/Astroide/Code/bin'

# Path to the working directory
WORKPATH='/Users/lacquema/Simulations/twa7/twa7_cb_a_dyn/twa7_cb_a_dyn_46'

# Generate the input files
$BINPATH/gen_multi_rmvs3 <<!
$BINPATH
rmvs3							# Method (rmvs3 / rmvs3_par)
0 10000000 0.05					# Initial, final, step time [yr]
10000 100000					# Output, dump frequencies [yr]
F T F F T F						# Options: int*2, real*4, nrj, debris jacobi coord, removing limits, J2 and J4 terms (T / F)
0.005 200. 10. -1. F			# Options: rmin, rmax, rmaxu, qmin, lclose (-1=not checked)
$WORKPATH
8								# Number of cpu cores per sub-simulation
1								# Units (0=[Msun,AU] / 1=[AU,yr])
1								# Coordinate (0=ecliptic / 1=invariant plane)
0								# Planet radius (0=None / 1=Physical / 2=Hill)
0.47							# Mass of central body 0 [Msun]
2								# Number of other bodies
0.34							# Mass of the orbiting body 1 [Mjup]
52.0 0.0 0.0 0.0 0.0 0.0	 	# Initial orbit parameters of 1: a[AU], e, i[deg], w[deg], W[deg], M[deg]
0.7								# Mass of the orbiting body 2 [Mjup]
17.0 0.0 0.0 0.0 0.0 0.0		# Initial orbit parameters of 2: a[AU], e, i[deg], w[deg], W[deg], M[deg]
654876543						# Random seed
200000							# Number of debris in the disk
25000							# Number of debris per sub-simulation
0. 0.01							# Acceptable range of eccentricity
2.								# Maximun of inclination
10. 40.							# Acceptable range of semi-major axis
200000							# Number of debris in the disk
25000							# Number of debris per sub-simulation
0. 0.01							# Acceptable range of eccentricity
2.								# Maximun of inclination
10. 40.							# Acceptable range of semi-major axis
0								# End of the inputs
!