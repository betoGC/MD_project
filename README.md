Jose A. Guerrero-Cruz\
fall2023-spring2024\
j.a.g.cruz@kjemi.uio.no\
albertogc1993@gmail.com
  
The script is based on the course-notes and diapositives kindly shared by Prof. Michelle-Cascella on multiscale modeling.

This is an experimental code I made to run locally in my computer, it is not ready for production calcs on the 
HPC systems, but you might want to play with it in your local-machine.

Generalities:\
Python script generated with version 3.8.2\
usage:
  
```bash
python3 MD.py
```
input files expected:

```python
  filename_pot='LC-10G-pot.data'     #Contains a parametrized potential in LC-gaussians read below
  filename_geo='fcc_ar.xyz'          #Contains a geometry, as an example a face-centered cubic lattice of solid argon (12 neighbors)
```
  
You can replace them with your method selection and geometries accordingly.
So-far we can handle only Argon structures, extension to water molecules soon.

Purpose is to perform a classical-MD simulation with a parameterized potential.
In this stage is possible only to work in cartesian geometries, transformations 
to internal coordinates or Z-matrix definition of the geometry could be benefitial 
to add bond-angle and dihedral-angle potentials. 

# Generating a parametrization potential for a pair AB

The script is capable of reading a linear combination of gaussians representing "pairs-of-atoms". 
You can generate such a description as follows: 

1. Swap a distance between particle A and B (assumed radial dependence only) to compute the total energy with a selected electronic-structure method. 
2. Then, set as zero the last point and shift the curve accordingly,
   
   $U_{AB}(r) = U_{r_i}-U_{r_n}$\
   such that\
   $U=0,\ F=0 \quad \forall \quad  r_{AB} > r_{n}$\
   
   this conditions impose a well in the local-minimum region and set to zero the force at the "infinite" limit.
  
3. Do a least-squares fitting of the numerical potential for a linear combination of gaussians with tempered exponents. (For sake of simplicity fit only the coefficients. The final expression of the potential will look similar to something like

   $U_{AB}(r) = \sum_{k} c_{k} e^{-z_{k}r_{AB}^2}$
   
   Non-linear fittings might be performed with the exponents $z_k$ , but I don't have a routine for that yet, make your choice.

4. Write the resulting coefficients and exponents into an ascii-file as if they were a XY coordinates in a molden format. It is not necessary to add an integer on the firstline with the number of gaussians followed by the coefficients and exponents.
   For example, having two gaussians, just add the coefficients and exponents:
   ```txt
     1.0  1.0
     0.5  0.3
   ```
   $kcal/mol$ units are assumed in the parameterized input potential, the energy is transformed to consistent units internally:
   
   $\left(\frac{8.314}{1.987}\right)\left(10^{-3}\right)\left(10^{-4}\right)$
   
   this is equivalent to the units transformation:
   
   $\frac{kcal}{mol} \rightarrow \frac{J}{mol} \rightarrow \frac{kg-Å^2}{ps^2-mol}$
   
6. Call The Force() function for pairwise forces calculated with the input potential using the 'gauss' flag.

# MD algorithm (limitations)

The script can handle a limited ammount of particles and still needs an efficient 
computational version. I'm pushing myself to vectorize the looping operations, probably its 
good idea switch to fortran or c++ and parallelized the code.
(I'm learning how to do it :D). 
  
Some arbitrary choices were made to reach a production-stage. Here some limitations.

## Algorithms for time integration:

  - Verlet
  - Leap-Frog
  - velocity Verlet

## Potentials for atom-pair forces 

  - Linear combination of gaussians $U_{AB} = \sum_k c_ke^{-z_kr_{AB}^2}$
  - Lennard-Jones $U_{AB} = \epsilon\left[\left(\frac{\sigma}{r_{AB}}\right)^{12}-\left(\frac{\sigma}{r_{AB}}\right)^{6}\right]$
  - Simple Harmonic-Oscillator $U= k(r_{AB}-r_{o})^2$  (Unreallistic k, and faulty PBCs don't use it yet for production calculations!)

## PBC in the minimal image convention

You need to define the Length of a cartesian centered squared box for PBC

The cutoff condition is defined as 
```python
if xx & yy & zz < 0.5*L:
  compute the Force...
else:
  Force=0.0
```
This assumes everything outside the box of length L as part of a non-interacting periodic image. Numerical instabilities near the edges were found with two and three-particles lying on a definition over let's say the z-axis with positions -3.0, 0.0, and 3.0 Å with a boxLenght of 6.0 Å. To circumvent the problem do: The net Force experienced by the ill-defined or outside the previous condition atom was set to 0.0. So, an artifitial motionless edge can be observed. After including random initial velocities compensated the instability problem. Crashing events were uncommon after this choice, remaining instability occur getting particles near the same position in the periodic image.

For our fcc-solid model of Argon, the box length include the farest atom and the first periodic image around 21.04 Å. This is computed after summation of the largest distance from the origin in one-axis and the first repeated atom.

## Thermostating

- Baerendsen
- Bussi-Donadio-Parinello
  
  Special attention is required to the $\ tau$, after experimenting with values between 0-5 ps, I found stability with a
  definition of 2.0 ps this is probably related to the time-steps selected around 5.0 fs which are three orders of magnitud
  different. However, the ratio between internal temperature and the set-point as well as the Wenier noise seem to compensate
  the numerical instability given square roots of numbers smaller than one. I need to look at the behavior for larger
  parameters. In any case, values larger than 1.0 ps cause the internal temperature to explode.

- A formal definition of a Wiener process is necessary for asignment of the stochastic term in the thermostatBDP function.
  The actual implementation scales a gaussian distribution for $N_f = 3*natom$ degrees of freedom between -0.1 to 0.1 centered
  around 0. In such a way, random sample velocities of a slow virtual distribution define a random walk similar to a Brownian
  process. Then, the collective variable to represent the Wiener noise filter was computed as:

  $W = \frac{1}{\sqrt{Nf}} \sum_k v_k \ \text{(random scaled velocity)}$
  
- About the initial conditions, random asignment of velocities with a 0 centered gaussian distribution with -1.0 to 1.0 spread-out (Maxwell-Boltzmann)

  $v(t_{0}) = w(t_{0})exp\left[ \frac{k_{b}T}{m} \right]$

  Removal of translational degrees of freedom for the center of mass as suggested in the notes is performed previously to every
  printing step to keep the total momentum nearby 0.0 even with thermostating on. This allows to keep the momentum below the
  third significant place in units of kg-Å/ps-gmol.
- A final note on the units, I'm thinking to switch into reduced variables to have a more friendly set
  The actual implementation have a very technical description in kg/gmol for mass, Å, ps, K and derived units for energy
  kg-Å**2/ps**2-K-gmol, all the proper transformations are computed in the corresponding places at the Forces calculations and
  before printing operations, all parameters and variables are consistent with this selection.

#TO-DO:
  1. Add a library list of elements and particles with properties like ionic-radius, gyration-Radii, VdW-Radius, mass, 1st, and 2nd ionization-Potentials,etc.
  2. Add the potentials already parameterized for the water-dimer. 
      This are nonBOA calculations at the Hartree-Fock APMO-MP2/aug-def2tzvp level calculations with a Nakai parametrized basis-set for protons and mono-dimensional in oxygen distances.
      Experimental Coarse-Grained motion of a water cluster is expected from such a potential description. The plan is try to compute energy of solvation for simple solutes.
  3. Compute the second virial coefficient from the classical-MD simulation of only waters using the inversion of the Boltzmann distribution.
      for a given temperature it must be an integral of the equilibrated energy profile.
  4. Think of computational efficiency, vectorization and parallelization.

# Disclaimer: 
This is an experimental code for a proof-of-principle idea: I want to compute a "coarse-grained" classical MD simulation for water to compute solvation energies, thermodynamic relations, and transport properties of water at room temperature and engineering operational regimes using PVT relations. I'm interested to plug-in a nonBOA potential developed by Flores-Moreno, et.al. around 2014 that was capable to handle hydrogens as "quantum" particles at the HF-MP2 level. In such an implementation, the PES of two-water molecules was reduced from a complex dimensional structure accounting for the two-water molecules orientation to a single-dimensional radial potential in the Oxygen-Oxygen distances because hydrogens were treated under the any-particle molecular orbital (APMO) simmilarly to handling of electrons in a restricted Hartree-Fock calculation assuming the hydrogens are behaving as fermions.

Despite some artifacts in the monodimensional PES, the parameterized potential allowed us to compute second-virial coefficients of water to some extent "reasonable", observing a better approximation including larger basis sets even if the stimated parameters were faulty in comparison to experimental data. Improving the potential description should lead to improve the accuracy of the method as well. With the goal to extend the applicability of these uncommon model for water to explicit solvation and dynamics here is my first attempt to connect both concepts: A simple model integrating classical-MD simulations with a non-BO-APMO HF-MP2 parameterized potential. 
