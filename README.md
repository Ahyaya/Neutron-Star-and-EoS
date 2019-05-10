# Neutron Star and Equation of State

### These codes are dedicated to compute non-rotating neutron stars with logarithmic equation of states. 

The available outputs in each project may include:

+ gravitational mass M
+ radius R
+ dimensionless moment of inertia  I
+ dimensionless tidal deformability  Lambda
+ baryon mass  Ma
+ pure mass  Mp

The current version of all projects use both testR-3 and testR-2 method to iterate instead of Runge-Kutta methods.

>## AScanRM
>This project uses 2-columns log(SI units) equation of states to compute R, M, I and Lambda. It's the very first code used in >our team which is barely able to automatically scan the central density then sort out the R-M-Lambda-I profile.
>
>We provide two modes to work:
>
>[0] The only configuration in auto-scan mode is the 2-columns log(SI units) equation of states, they must be correctly placed in the directory ~/AScanRM/EoS_lib in the form of \*.txt files.
>
>[1] In addition to EoS files, you must determine the central density sequence in ~/AScanRM/Central_Density/DenSeq.txt in the units of kg/m^3

## Ad_GW170817
Likelihood estimation from GW170817 using the built-in ANM EoS.

The two versions (LE_v10 and LE_Mq) are slightly different in their likelihood function, but overall, they provide posterior distribution by matching the ANM parameters to GW170817.

The \*.m files are for octave-4.4.1 or higher version.
