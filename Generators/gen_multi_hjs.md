# Documentation to adapt the `gen_multi_hjs.sh` script

This script generates sub-simulations input files for the `HJS` integrator, which is designed to simulate the dynamics of planetary systems containing debris disks and bodies of comparable mass. The generated files specify all necessary simulation parameters in the required order. Below is a description of each parameter, listed in the order they appear in the input file:

**Path to executables**
   - *Format*: `str`
   - *Description*: Absolute path indicates where to find the required executables.

**Integration method**
   - *Value*: `hjs` or `hjs_par`
   - *Description*: Choice of integration method (serial or parallel).

**Simulation times**
   - *Format*: `t_0` `t_f` `dt`
   - *Type*: `float` `float` `float`
   - *Description*: 
        - `t_0`: Start time of the simulation (in years).
        - `t_f`: End time of the simulation (in years).
        - `dt`: Integration time step (in years). Choose it carefully to balance accuracy and performance; for example, `0.05` means the step is 5% of the shortest orbital period in the system.

**Output frequencies**
   - *Format*: `output_freq` `dump_freq`
   - *Type*: `float` `float`
   - *Description*: 
      - `output_freq`: Output frequencie (in years).
      - `dump_freq`: Dump frequencie (in years).

**Computational options**
   - *Format*: `F1` `F2` `F3` `F4` `F5` `F6`
   - *Type*: `Bool` `Bool` `Bool` `Bool` `Bool` `Bool`
   - *Description*: If `True`,
       - `F1`: Use 2-byte integers (saves memory, may limit range),
       - `F2`: Use 4-byte floating-point numbers (reduces memory usage, may affect precision).
       - `F3`: Compute and output system energy.
       - `F4`: Use Jacobi coordinates for debris particles.
       - `F5`: Enable removal limits (activates the removal criteria section below).
       - `F6`: Include J2 and J4 gravitational moments (for oblate central bodies).

**Removal limits** (if removal limits are enabled)

- *Format*: `rmin` `rmax` `rmaxu` `qmin` `lclose`
- *Type*: `float` `float` `float` `float` `bool`
- *Description*: (All distances in AU; use `-1` to disable a specific check)
    - `rmin`: Remove particles if they come closer than this minimum distance from the central body.
    - `rmax`: Remove particles if they move farther than this maximum distance from the central body.
    - `rmaxu`: Upper limit for the central body's radius when not gravitationally bound.
    - `qmin`: Remove particles if their perihelion distance falls below this value.
    - `lclose`: If `True`, remove particles upon close encounter with a planet (encounter radius is specified in the `io_init_pl` file).

**Working directory**
   - *Type*: `str`
   - *Description*: Absolute path indicates where running files will be generated.

**Number of CPU cores per sub-simulation**
   - *Type*: `int`
   - *Description*: Number of CPU cores used for each sub-simulation.

**Units**
   - *Values*: 
      - `0`: Mass in Msun, distance in AU.
      - `1`: distance in AU, time in years.
   - *Description*: Sets the units used for all physical quantities in the simulation.

**Coordinate system**
   - *Values*: 
      - `0`: Ecliptic plane.
      - `1`: Invariant plane.
   - *Description*: Determines the reference plane for all orbital elements. The ecliptic plane is the plane of Earth's orbit, while the invariant plane is defined by the total angular momentum of the system.

**Planetary radius type**
   - *Values*: 
      - `0`: No planetary radius considered.
      - `1`: Use the physical radius of the planet.
      - `2`: Use the Hill radius of the planet.
   - *Description*: Specifies how the radius of each planet is defined for removal criteria and close encounter calculations.

**Number of bodies**
   - *Type*: `int`
   - *Description*: Number of planets or massive objects orbiting the central body.

**Bodies mass**

For each body, repeat the next line:

1. **Mass of the body** 
   - *Type*: `float`
   - *Description*: Mass of the body in Msun. The first will be the central mass. 

**Orbiting bodies definitions**

For each orbiting body, repeat the next 2 lines:

1. **Orbital hierarchy of the orbiting body**
    - *Format*: `B1` ... `Bn` (where `n` is the number of bodies)
    - *Values*: For each body, specify:
        - `1`: The body is the reference of the hierarchy.
        - `-1`: The body is inside the orbit of the reference body.
        - `0`: The body is outside the orbit of the reference body.
    - *Description*: Defines the hierarchical structure of the system by specifying, for each body, its relationship to the others. 

2. **Initial orbital parameters of the orbiting body**
    - *Format*: `a` `e` `i` `w` `W` `M`
    - *Type*: `float` `float` `float` `float` `float` `float`
    - *Description*:
      - `a`: Semi-major axis in AU.
      - `e`: Eccentricity.
      - `i`: Inclination in deg.
      - `w`: Argument of perihelion in deg.
      - `W`: Longitude of ascending node in deg.
      - `M`: Mean anomaly in deg.

**Random seed**
   - *Type*: `int`
   - *Description*: Seed for initializing the random number generator.

**Debris disk definitions**

For each debris disk, repeat the next 7 lines:

1. **Number of debris in the disk**
   - *Type*: `int`
   - *Description*: Total number of test particules in the disk.

2. **Number of debris per sub-simulation**
   - *Type*: `int`
   - *Description*: Number of test particles assigned to each sub-simulation. This parameter controls the parallelization: the total debris disk will be split into multiple sub-simulations, each containing this number of particles.

3. **Orbital hierarchy of the disk**
   - *Format*: `B1` ... `Bn` (where `n` is the number of bodies)
   - *Values*: For each body, specify:
      - `-1`: The body is inside the disk.
      - `0`: The body is outside the disk.
   - *Description*: Defines the hierarchical structure of the system by specifying, for each disk, its relationship to the bodies. 

4. **Reference frame of the disk**
    - *Values*: 
        - `-1`: Aligned with the orbital plane of the reference body.
        - `0`: Aligned with the ecliptic plane.
        - `1`: Aligned with the system's center of mass.
        - `2`: Aligned with the invariant plane.
    - *Description*: Specifies the coordinate system in which the disk is initialized.

5. **Eccentricity range in the disk**
   - *Format*: `e_min` `e_max`
   - *Type*: `float` `float`
   - *Description*: 
      - `e_min`: Minimum eccentricity for particles in the disk.
      - `e_max`: Maximum eccentricity for particles in the disk.

6. **Maximum inclination in the disk**
   - *Type*: `float`
   - *Description*: Maximum inclination for particules in the disk in deg.

7. **Semi-major axis range in the disk**
    - *Format*: `a_min` `a_max`
    - *Type*: `float` `float`
    - *Description*: 
        - `a_min`: Minimum semi-major axis for particles in the disk in AU.
        - `a_max`: Maximum semi-major axis for particles in the disk in AU.

**End of input**
   - *Value*: `0`
   - *Description*: Indicates the end of the input parameters.


