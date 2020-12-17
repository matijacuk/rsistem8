# rsistem8

This is the numerical integrator used in the manuscript by Cuk, Lock, Stewart and Hamilton submitted to PNAS in Dec 2020. The FORTRAN file is the same in all simulations, and the input files are for the state 1 Myr into the simulation shown in Fig 1, which is when we switch to the standard timestep and before the Laplace Plane instability starts. The basics of numerical methods are described in the paper's SI file, and the parts that are the same as in Cuk et al. 2016 Nature paper are described in the SI of that paper. Here we will explain the format of the input files. Fortmat of someof them is shared with the previous integrator SIMPL used for dynamics of giant planet satellite systems.

Units used throughout are AU and years. 

At recording intervals, the program makes .rec files which have the same format as .in files and can be used for restrart

settings.in:
Line 1: beginning time, end time, timestep
Line 2: Recording interval, output interval
Line 3: five flags [0=NO, 1=YES]: 
1. append .out files?
2. include direct planetary perturbations on the satellites [irrelevant in this paper]
3. include tidal forces?
4. include artificial migration [not used in RSISTEM8, keep 0]
5. output format [0=using Z=0 reference plane and AU, 1=using planetary radii and equator, 2=uising planetary equator and heliocentric orbit]

planets.in:
Line 1: Solar mass (=4 pi^2), N of planets, input format flag [0=heliocentric vectors, 1=heliocentric elements] 
Line 2: Earth's mass relative to the Sun
Line 3: x, y, z heliocenytic positions of Earth in AU
line 4: vx, vy, vz helocentric velocities on Earth in AU/year

moons.in
Line 1: number of the planet with moons, number of moons, input format flag [0=planetocentric vectors, 1=planetocentric elements]
Line 2: (J_2 R^2) for Earth, x, y, z components of the spin axis unit vector
Line 3: lunar mass relative to solar
Line 4: x, y, z planetocentric positions of the Moon
line 5: vx, vy, vz planetocentric 

tidal.in:
Line 1: Earth's radius, Earth's tidal Q, Earth's tidal Love number, Earth's spin rate [in 2pi/years]
Line 2: Lunar radius, lunar tidal Q, lunar tidal Love number

precess.in: only one value, dimensionless moment of inertia of Earth at low spins. 

lunar.in:
Line 1: number of the satellite with resolved rotation
Line 2: 3 components of lunar rotation around 3 principal axes
Line 3: Lunar dimensionless moments of inertia around 3 principal axes
Lines 4-6: 9 elements of the conversion matrix between lunar body frame and the inertial frame

The output creates (or appends) the following output files:

planet101: time, Earth's heliocentric elements (a in AU, e, inc, long. node, arg. perihelion, mean anomaly)

moon101: time, lunar planetocentric orbital elements (a, e, inc, long. node, arg. perigee, mean anomaly) [units depend on the last flag in settings.in]

pole.out: time, Earth's obliquity wrt. ecliptic, long. node, spin rate

spin_orb.out: time, lunar obliquity wrt geocentric orbit, associated long. node, spin rate

spin_ecl.out: time, lunar obliquity wrt Earth's heliocentric orbit, associated long. node, spin rate

long.out: time, 3 spin rates around lunar principal axes, angle between the lunar longes axis and Earth-Moon line [in radians]











