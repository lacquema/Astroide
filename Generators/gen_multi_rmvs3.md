# Parameter Documentation for the gen_multi_hjs.sh Script

This script generates the input files for the HJS integrator. Each parameter is described below, based on the comments present in the script.

- **BINPATH**: Path to the required executables (e.g., gen_multi_hjs, hjs, etc.).
- **WORKPATH**: Path to the working directory where the input files will be generated.

The parameters passed to the gen_multi_hjs executable are:

1. **BINPATH**  
    Path to the executables.

2. **Method**  
    `hjs` or `hjs_par`: integration method to use.

3. **Initial time, final time, time step [yr]**  
    Format: `0. 10000000. 0.05`  
    Start time, end time, and time step in years.

4. **Output and dump frequencies [yr]**  
    Format: `10000. 100000.`  
    Output frequency and dump (save) frequency in years.

5. **Options**  
    Format: `F T F F T F`  
    - int*2: Use 2-byte integers (T/F)
    - real*4: Use 4-byte reals (T/F)
    - nrj: Compute energy (T/F)
    - debris jacobi coord: Jacobi coordinates for debris (T/F)
    - removing limits: Particle removal limits (T/F)
    - J2 and J4 terms: Include J2 and J4 terms (T/F)

6. **Removal limits (if enabled)**  
    Format: `0.005 200. 10. -1. F`  
    - rmin: Minimum radius
    - rmax: Maximum radius
    - rmaxu: Ultimate maximum radius
    - qmin: Minimum perihelion (-1 = not checked)
    - lclose: Enable close removal (T/F)

7. **WORKPATH**  
    Path to the working directory.

8. **Number of CPU cores per sub-simulation**  
    Integer value.

9. **Units**  
    0 = [Msun, AU], 1 = [AU, yr]

10. **Coordinate plane**  
     0 = ecliptic, 1 = invariant plane

11. **Planetary radius**  
     0 = None, 1 = Physical, 2 = Hill radius

12. **Number of bodies**  
     Integer value.

13. **Mass of body 1 [Msun]**  
     Floating-point value.

14. **Mass of body 2 [Msun]**  
     Floating-point value.

15. **Mass of body 3 [Msun]**  
     Floating-point value.

16. **Orbital hierarchy of body 2**  
     Format: `-1 1 0`  
     -1 = inside, 0 = outside, 1 = body in question

17. **Initial orbital parameters of body 2**  
     Format: a[AU], e, i[deg], w[deg], W[deg], M[deg]  
     - a: semi-major axis [AU]
     - e: eccentricity
     - i: inclination [deg]
     - w: argument of perihelion [deg]
     - W: longitude of ascending node [deg]
     - M: mean anomaly [deg]

18. **Orbital hierarchy of body 3**  
     Same format as for body 2.

19. **Initial orbital parameters of body 3**  
     Same format as for body 2.

20. **Random seed**  
     Integer used to initialize the random number generator.

21. **Number of debris in disk 1**  
     Integer value.

22. **Number of debris per sub-simulation (disk 1)**  
     Integer value.

23. **Orbital hierarchy of disk 1**  
     Format: `-1 0 0`  
     -1 = inside, 0 = outside

24. **Coordinate plane of disk 1**  
     -1 = orbit, 0 = ecliptic, 1 = centers, 2 = invariant plane

25. **Acceptable eccentricity range (disk 1)**  
     Format: `0. 0.01`

26. **Maximum inclination (disk 1)**  
     Floating-point value.

27. **Acceptable semi-major axis range (disk 1)**  
     Format: `5. 15.`

28. **Number of debris in disk 2**  
     Integer value.

29. **Number of debris per sub-simulation (disk 2)**  
     Integer value.

30. **Orbital hierarchy of disk 2**  
     Same format as for disk 1.

31. **Coordinate plane of disk 2**  
     Same format as for disk 1.

32. **Acceptable eccentricity range (disk 2)**  
     Format: `0. 0.01`

33. **Maximum inclination (disk 2)**  
     Floating-point value.

34. **Acceptable semi-major axis range (disk 2)**  
     Format: `15. 130.`

35. **End of inputs**  
     `0`: marks the end of the input parameters.