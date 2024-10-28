import matplotlib.pyplot as plt
import csv
import numpy as np
import glob
import random
import sys

#Stablishing object for random normal distributed numbers generator
rng = np.random.default_rng()

#MD algorithm. Initialization of simulation variables starts here
#MD_time [ps] in general the time units used in the simulation are assumed as pico-seconds
#make all the relevant unit conversions consistent with: kg,Å,ps,kg/mol,Å/ps,etc.
#kb units assumed in kg-Å**2/ps**2-K-mol
#dt [fs]
#nsteps is an integer ratio
MD_time = 20          ### ps    
dt      = 5           ### fs    
dt = dt/1000.0        # fs -> ps        
nsteps  = MD_time/dt      
nsteps  = int(nsteps)

###Size of the cartesian Box for periodic-boundary-conditions in the minimal image formalism ###
L  = 21.04  ##Å
#L = 9.0    ##Å
print("MD_steps",nsteps,"Tot_Sim_Time [ps] =%10.6f"%(nsteps*dt),"Time_Step [fs]=%10.6f"%(dt*1000))
print("PBC boxLenght:%10.6f [Å]"%(L))

#flag: kind of pair-potential, Vij and Fij: 
#LJ:     6-12 Lennard-Jones,
#harm:   harmonic-oscilator,
#gauss:  Numerical potential fitted to a linear combination of gaussians
#flagti: Time integration algorithm: 
#Verlet: Classical Verlet algorithm for integration of time
#leapFrog: Time integration adding acceleration from previous step 
#velVerlet: Verlet modified to include acceleration and forces each step with randomized initial distributions.

#Recommended
#flag='gauss'
#flagti='velVerlet'

flag='LJ'
flagti='velVerlet'

#Input files, 
#filename_pot: A parametrized potential as a LC of gaussian functions
#filename_geo: Cartesian xyz geometry in molden std format

if flag == 'harm':
  filename_pot='LC-10G-pot.data'
  filename_geo='fcc_ar.xyz'
  #filename_geo='test.xyz'
elif flag == 'LJ':
  filename_pot='LC-10G-pot.data'
  #filename_geo='LJ-1.xyz'
  filename_geo='fcc_ar.xyz'
elif flag == 'gauss':
  #filename_pot='Gauss-ccpvtz.dat'
  filename_pot='Gauss-ccpvdz.dat'
  filename_geo='testO-15Å.xyz'

#Reading geometry from a molden .xyz format
#cartesian coords in Å 
#N: number of atoms 
#coords:[iatom,[x,y,z]]
def readGeom(filexyz):
  f = filexyz
  xyz = open(f)
  N = int(xyz.readline())
  header = xyz.readline()
  atom_symbol, coords = ([] for i in range (2))
  for line in xyz:
      atom,x,y,z = line.split()
      atom_symbol.append(atom)
      coords.append([float(x),float(y),float(z)])
  xyz.close()
  return N,atom_symbol,coords

#Reading the c_k and z_k for the LC-gaussians potential
#Uij=\sum_k c_k*exp(-z_k*rij**2)
def readPot(filename):
  f2 = open(filename, 'r')
  lines = f2.readlines()
  x = []
  y = []
  for line in lines:
      p = line.split()
      x.append(float(p[0]))
      y.append(float(p[1]))
  f2.close()
  return x,y

#Calculation of the radial gradient given a LC-gaussians as a potential
#for a pair of particles.# Analytical derivative given a one-dimensional gaussian
#Minimal image convention on a cartesian cube
def Force(atoms,coords,c_k,z_k,L):
  nGauss=len(c_k)
  pos1=np.zeros(3)
  pos2=np.zeros(3)
  F=np.zeros((natom,3))
  U=np.zeros(natom)
  F[:,:] = 0.0
  Fr=0.0
  for i in range(natom-1): #Loop over atom i. 
    if atoms[i] == 'H':
      continue
    else:
      pos1[0]=coords[i,0]
      pos1[1]=coords[i,1]
      pos1[2]=coords[i,2]
      for j in range(i+1,natom): # Loop starts at i+1 atom, and finish in natom-1. For two atoms this must be only the 1 position.
        if atoms[i] == 'H':
          continue
        else:
          pos2[0]=coords[j,0]
          pos2[1]=coords[j,1]
          pos2[2]=coords[j,2]
          xx  = pos2[0]-pos1[0]           # distance between 1 and 2 on cartesian component i
          yy  = pos2[1]-pos1[1]
          zz  = pos2[2]-pos1[2]
          nLx = np.round(xx/L)
          nLy = np.round(yy/L)
          nLz = np.round(zz/L) 
          xx -= L*nLx
          yy -= L*nLy
          zz -= L*nLz
          rr = xx**2+yy**2+zz**2 # Squared distance from the image to the iatom 
          r  = np.sqrt(rr)
          if (xx < 0.5*L) & (yy < 0.5*L) & (zz < 0.5*L): #Cutoff condition, half the cell length
            Fr   = r*np.dot(c_k[:]*z_k[:],np.exp(-rr*z_k))
            Fr  *= -2.0*8.314e-4*1.0e-3/1.987
            acum = np.sum(Fr)
            U[i] = np.sum(c_k[:]*np.exp(-rr*z_k))
            F[i,0] += -acum*xx/(r)
            F[i,1] += -acum*yy/(r)
            F[i,2] += -acum*zz/(r)
            F[j,0] -= F[i,0]      #Mirroring the force to the j-th particle once i-th computed
            F[j,1] -= F[i,1]
            F[j,2] -= F[i,2]
          else:
            F[i,0] += 0.0
            F[i,1] += 0.0
            F[i,2] += 0.0
            F[j,0] -= F[i,0]      #Mirroring the force to the j-th particle once i-th computed
            F[j,1] -= F[i,1]
            F[j,2] -= F[i,2]
            continue
      #U[iatom] += Ur  
      #print("iatom,jatom",iatom,jatom,"R",r,"Fi=",F[iatom],"Fj=",F[jatom])
  #print("Force array=\n",F)
  #sys.exit("code debug stop")
  return F,coords,U

#Calculating the force for a harmonic-oscilator potential
def ForceH(coords,natom,L):
  xe = 1.5
  k  = 16.0         #J/mol
  pos1=np.zeros(3)
  pos2=np.zeros(3)
  F=np.zeros((natom,3))
  U = np.zeros(natom)
  F[:,0:3] = 0.0
  for iatom in range(natom-1):
    pos1=coords[iatom,0:3]
    for jatom in range(iatom+1,natom):
      pos2=coords[jatom,0:3]
      xx = pos2[0]-pos1[0]
      yy = pos2[1]-pos1[1]
      zz = pos2[2]-pos1[2]
      nLx = np.round(xx/L)
      nLy = np.round(yy/L)
      nLz = np.round(zz/L) 
      xx = xx-L*nLx
      yy = yy-L*nLy
      zz = zz-L*nLz
      r=np.sqrt(xx**2+yy**2+zz**2)
      U[iatom]   += 0.5*k*(r-xe)**2
      F[iatom,0] += -k*(r-xe)*xx/r
      F[iatom,1] += -k*(r-xe)*yy/r
      F[iatom,2] += -k*(r-xe)*zz/r
      F[jatom,0]+= -F[iatom,0]
      F[jatom,1]+= -F[iatom,1]
      F[jatom,2]+= -F[iatom,2]
  return F

def ForceLJ(coords,natom,L,T):
  #For a given molecular geometry compute the force(natom) array
  N  = 6.022e23
  kb = 8.314e-4             # kg-Å2/ps2-K-mol
  epsilon_over_k = 122.4   # K^-1  From R. Byron Bird, Transport phenomena, pp. 277 
  sigma = 3.432            # Å
  epsilon = epsilon_over_k*kb
  F=np.zeros((natom,3))
  U=np.zeros(natom)
  pos1=np.zeros(3)
  pos2=np.zeros(3)
  Fr=0.0
  for i in range(natom-1):        # Visits only the upper triangular Hessian and mirror the pairs force
    pos1[0] = coords[i,0]         # Rows going to the last-1 element
    pos1[1] = coords[i,1]
    pos1[2] = coords[i,2]
    for j in range(i+1,natom):    # Columns going from the diagonal to the matrix end
      pos2[0] = coords[j,0]
      pos2[1] = coords[j,1]
      pos2[2] = coords[j,2]
      xx  = pos2[0]-pos1[0]
      yy  = pos2[1]-pos1[1]
      zz  = pos2[2]-pos1[2]
      nLx = np.round(xx/L)
      nLy = np.round(yy/L)
      nLz = np.round(zz/L) 
      xx -= L*nLx
      yy -= L*nLy
      zz -= L*nLz
      rr  = xx**2+yy**2+zz**2
      ###add cutoff restriction
      if (xx < 0.5*L) & (yy < 0.5*L) & (zz < 0.5*L):
        X   = 12.0*(sigma**2/rr)**6
        Y   = 6.0*(sigma**2/rr)**3 
        U[i] += 4.0*epsilon*((sigma**2/rr)**6-(sigma**2/rr)**3)
        Fr  = 4.0*epsilon*(X-Y)    # Pairwaise radial gradient 
        F[i,0] += -Fr*xx/rr        # Adding in the i-th atom
        F[i,1] += -Fr*yy/rr
        F[i,2] += -Fr*zz/rr
        F[j,0] -= F[i,0]           # Mirroring the j-th contribution
        F[j,1] -= F[i,1]
        F[j,2] -= F[i,2]
      else:                        #Warning: Forces contributions set to zero for atom pairs outside the cutoff condition
        F[i,0] += 0.0   
        F[i,1] += 0.0
        F[i,2] += 0.0
        F[j,0] -= F[i,0]           # Mirroring the j-th zero contribution
        F[j,1] -= F[i,1]
        F[j,2] -= F[i,2]
        continue
      ###########################
      #print("iatom",i,"jatom",j,"xx %10.6f,yy %10.6f, zz%10.6f"%(xx,yy,zz),"rr=%10.6f"%(rr),"[nLx,nLy,nLz]=",nLx,nLy,nLz,"Fr=%10.6f"%(Fr))
  #sys.exit("code debug stop")
  return F,U

def removeTranslation(natom,m,v):
  p_com = np.zeros(3)
  Mtot =np.sum(m)
  v[:,0] -= np.sum(m[:]*v[:,0])/Mtot
  v[:,1] -= np.sum(m[:]*v[:,1])/Mtot
  v[:,2] -= np.sum(m[:]*v[:,2])/Mtot
  
  return v
    
def thermostatB(natom,dt,v,T_o,int_T):
    tau = 2.0 #ps
    B=(dt/tau)*(1.0-T_o/int_T)
    lambda_B = np.sqrt(1.0-B)
    new_v=np.zeros((natom,3))
    new_v = lambda_B*v
    #print("lambda_B=",lambda_B)
    #print("dt",dt,"tau",tau,"T_o",T_o,"int_T",int_T)
    #print("old_v",v)
    #print("new_v",new_v)
    #sys.exit("code debug stop")
    return new_v

def thermostatBDP(natom,dt,v,T_o,int_T,m):
  tau = 2.0 #ps
  Nf=3*natom
  kb=8.314e-4
  w=0.0
  for i in range(Nf):
    w+=np.sqrt(kb*T_o/m)*random.uniform(-0.10000,0.10000)
  w*=1.0/np.sqrt(Nf)
  B  = 1.0-(dt/tau)*(1.0-T_o/int_T)
  St = np.sqrt((4.0/3.0/natom)*(T_o/int_T))*(w/np.sqrt(tau))
  lambda_B = np.sqrt(B+St)
  v*=lambda_B
  return v

#Reading the external potential
z,c=readPot(filename_pot)
c=np.asarray(c)
z=np.asarray(z)
#print("Here is the coefficient",c)
#print("Here is the zeta",z)
#sys.exit("code debug stop")

#Initialization of coordinates 
natom,atoms,coords=readGeom(filename_geo)
coords=np.asarray(coords)

#Uncomment for validation of numbers 
#print(c,z)
#print("This is the initial geometry:")
#print(natom)
#for i in range(natom):
#  print(atoms[i],"%20.10f"%(coords[i,0]),"%20.10f"%(coords[i,1]),"%20.10f"%(coords[i,2]))

#Initialization of physical attributes 
m=np.zeros(natom)              #Mass
N    = 6.022e23                #Avogadro's constant - particles/mol missused
if flag == 'gauss':
  for i in range(natom):
    if atoms[i] == 'O':
      #m[i] = 15.9994e-3
      m[i] = 18.0152e-3
    elif atoms[i] == 'H':
      m[i] = 1.0079e-3  
    elif atoms[i] == 'Ar':
      m[i] = 39.948e-3     #Mass in - Kg/mol Assumed Water particles. 
elif flag == 'harm':           #TO-DO: Replace with a dictionary of particles 
  for i in range(natom):
    if atoms[i] == 'O':
      m[i] = 15.9994e-3
    elif atoms[i] == 'H':
      m[i] = 1.0079e-3  
    elif atoms[i] == 'Ar':
      m[i] = 39.948e-3     #Mass in - Kg/mol Assumed Water particles. 
elif flag == 'LJ':
  for i in range(natom):
    if atoms[i] == 'O':
      m[i] = 15.9994e-3
    elif atoms[i] == 'H':
      m[i] = 1.0079e-3  
    elif atoms[i] == 'Ar':
      m[i] = 39.948e-3     #Mass in - Kg/mol Assumed Water particles. 

################################Allocation of physical variables and parameters for units consistency########################
v=np.zeros((natom,3))          #Velocities
p=np.zeros((natom,3))          #Boltzmann probabilities on each cartesian component
freq=np.zeros(natom)           #Boltzmann Frequencies 
F=np.zeros((natom,3))          #Force array in kg-Å**2/ps**2 instead of kg-m**2/s**2 = Joules-m
U    =np.zeros(natom)          #Potential Energy, remember 1st law: the total internal energy E = K +/- U = Q +/- W be aware of consistency 
pi   = 2.0*np.arccos(0.0)      #pi constant intrinsic definition in numpy 
kb   = 8.314e-4                #Boltzmann constant Kg-Å**2/ps**2-K-mol  (handling of numbers/encourage to use reduced units)
T    = 25.0                    #Assumed absolute temperature of the system in Kelvin
Ekin = (3.0/2.0)*kb*T          #Initial kinetic energy accordingly to kinetic theory of gases, Kg-Å**2/ps***2-mol

print("Potential:",flag)
print("Time-Integration:",flagti)
print("Init_kin_Energy=%10.6f"%(Ekin),"[kg-Å**2/ps**2-mol]")
print("meanVelocity=%10.6f"%(np.sqrt(2.0*Ekin/m[0])),"[Å/ps]")
print("Thermostat SetPoint/Initialization Temperature=%10.6f"%(T),"[K]")
print("Trajectory-coordinates file: a.trj")
print("Forces, internal energy: F.trj")
print("velocities: vel.trj")
print("time-steps vs internal Temperature: a.txt")

###############################Initialization of velocities###########################################
for i in range(natom):
  for j in range(3):
    v[i,j]=np.sqrt(kb*T/m[i])*random.uniform(-1.000,1.000)
    #print("iatom",i,"arg",arg)

##############################Removing translational motion of the center of mass ###################
v=removeTranslation(natom,m,v)
#print("======Initial set of velocities======")
#print("atom \t vx \t vy \t vz")
#for i in range(natom):
#  print(atoms[i], "%10.6f %10.6f %10.6f"%(v[i,0],v[i,1],v[i,2]))
####################Finalization of setting up initial conditions #################################

wcoords=coords
newcoords=coords
#v[:,:] = 0.0
totalForce=np.zeros(3)
totalP=np.zeros(3)
av_V2=np.zeros(natom)
with open('a.trj', 'w') as trj:            # MD trajectory 
  with open('F.trj', 'w') as fbwu:         # Computed Forces by atom
    with open('vel.trj', 'w') as vel:      # Computed velocities by atom
      with open('a.txt', 'w') as report:   # report internal temperature and kinetic energy in xy format for gnuplot post-processing

        #MasterLoop in time#####
        for t in range(nsteps+1):

          #Calculate the force vector, remember each particle experience force due to all possible pairs 
          if flag == 'gauss':
            F,wcoords,U=Force(atoms,wcoords,c,z,L)
          elif flag == 'harm':
            F=ForceH(wcoords,natom,L)
          elif flag == 'LJ':
            F,U=ForceLJ(wcoords,natom,L,T)
         
          #Veryfying the total momentum and forces on each step, force is correct, momentum has variations on 10e-15 without thermostat 
          T_o=T         #thermostat set-point To = T to asign random initial velocities
          eta = 3.0     #dimensionality for thermostating we are in a 3D cartesian definition
          totalP[:] = 0.0
          
          av_V2  = v[:,0]*v[:,0]  
          av_V2 += v[:,1]*v[:,1]  
          av_V2 += v[:,2]*v[:,2]  
          av_K = np.dot(m,av_V2)
          #print(av_K)
          #sys.exit("code debug stop")
      
          totalForce[0] = np.sum(F[:,0])
          totalForce[1] = np.sum(F[:,1])
          totalForce[2] = np.sum(F[:,2])

          totalP[0] = np.dot(m,v[:,0])
          totalP[1] = np.dot(m,v[:,1])
          totalP[2] = np.dot(m,v[:,2])

          ##Internal temperature of the actual time step T(t), Kinetic energy K(t), potential energy U(t), and internal energy E(t)#
          int_T=av_K/(eta*natom*kb)
          Ekin=0.5*av_K
          Epot=np.sum(U)
          E=Ekin+Epot

          #Printing headers for trajectory files each 10 time-steps
          if t%10  == 0:
            print(natom, file=trj)
            print("nstep:",t,"time[ps]= %10.6f"%(dt*t),\
            "int_T [K]= %10.6f"%(int_T),\
            "Ekin [kJ/mol]= %10.6f"%(Ekin*1e1),\
            "Epot [kJ/mol]= %10.6f"%(Epot*1e1),\
            "E [kJ/mol]= %10.6f"%(E*1e1),file=trj)
            
            print(natom, file=fbwu)
            print("nstep:%5d"%t,"time[ps]= %10.3f"%(dt*t),\
                  "TotalForce: %10.6f %10.6f %10.6f"%(totalForce[0],totalForce[1],totalForce[2]),\
                  "Total momentum: %10.6f %10.6f %10.6f"%(totalP[0],totalP[1],totalP[2]),\
                  "Int_T[K]= %10.6f"%(int_T),file=fbwu)
            print(natom, file=vel)
            print("nstep:%5d"%t,"time[ps]=%10.3f"%(dt*t),file=vel)
            #print("%10.6f %10.6f %10.6f"%(dt*t,int_T,0.5*av_K),file=report)

            #Removing translation of mass center every t%x steps; TO-DO: add as input parameter nsteps for printing in the modulus conditions#######
            v=removeTranslation(natom,m,v)

          ######Time integration algorithms, try to improve the case structure into a single flag comparison########################

          for iatom in range(natom):
            if t%10 == 0:
              print(atoms[iatom]," ","%20.10f"%(wcoords[iatom,0])," ","%20.10f"%wcoords[iatom,1]," ","%20.10f"%wcoords[iatom,2], file=trj)
              print(atoms[iatom]," ","%10.6f"%(F[iatom,0])," ","%10.6f"%(F[iatom,1])," ","%10.6f"%(F[iatom,2]), file=fbwu)
              print(atoms[iatom]," ","%10.6f"%(v[iatom,0])," ","%10.6f"%(v[iatom,1])," ","%10.6f"%(v[iatom,2]), file=vel)
            for j in range(3):
               if flagti == 'Verlet':
                 if t == 0:
                   newcoords[iatom,j] = wcoords[iatom,j]+F[iatom,j]*dt**2/m[iatom]
                   oldcoords=wcoords
                   wcoords=newcoords 
                 else:
                   newcoords[iatom,j] = 2.0*wcoords[iatom,j]-oldcoords[iatom,j]+(F[iatom,j]/m[iatom])*dt**2
                   oldcoords=wcoords
                   wcoords=newcoords 
               elif flagti == 'leapFrog':
                 if t == 0:
                   newcoords[iatom,j] = wcoords[iatom,j] + v[iatom,j]*dt
                   v[iatom,j] = v[iatom,j] + F[iatom,j]*dt/m[iatom]
                   oldcoords=wcoords
                   wcoords=newcoords 
                 else:
                   #v=thermostatB(natom,dt,v,T_o,int_T)                #It works, TO-DO: Implement a flag selection of the thermostat in a single call with none, B, and BDP options
                   v=thermostatBDP(natom,dt,v,T_o,int_T,m[iatom])      #Warning, mass of the i-th atom; replace by mass of the i-th particle
                   newcoords[iatom,j] = wcoords[iatom,j] + v[iatom,j]*dt
                   v[iatom,j] = v[iatom,j] + F[iatom,j]*dt/m[iatom]
                   oldcoords=wcoords
                   wcoords=newcoords 
               elif flagti == 'velVerlet':
                 if t == 0:
                   oldForce=F
                   newcoords[iatom,j] = wcoords[iatom,j] + v[iatom,j]*dt + 0.5*F[iatom,j]*dt**2/m[iatom]
                   v[iatom,j] = v[iatom,j] + F[iatom,j]*dt/m[iatom]
                   oldcoords=wcoords
                   wcoords=newcoords 
                 else:
                   newcoords[iatom,j] = wcoords[iatom,j] + v[iatom,j]*dt + 0.5*F[iatom,j]*dt**2/m[iatom]
                   dF = F[iatom,j] + oldForce[iatom,j]
                   #v=thermostatB(natom,dt,v,T_o,int_T)                #It works
                   v=thermostatBDP(natom,dt,v,T_o,int_T,m[iatom])    #Warning, mass of the i-th atom/particle
                   v[iatom,j] = v[iatom,j] + 0.5*dF*dt/m[iatom]
                   oldcoords=wcoords
                   wcoords=newcoords 
        trj.close()
        fbwu.close()
        vel.close()
        report.close()
