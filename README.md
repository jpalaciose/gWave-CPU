# gWave-CPU
CPU-based parallel tsunami simulation code for a single domain.

# How to use
1. Install Fortran compiler with MPI (it can be downloaded from OneAPI website)
2. Compile the code with: mpif90 gWave-CPU_mp.f90 gWave-CPU_ut.f90 -O3 -o gWave-CPU.exe
3. For the simulation test, run with 48 parallel processes: mpirun -n 48 ./gWave-CPU.exe
4. The propagation results will be written in the directory: 2_results
   
Note that the program is defined for a single domain. In future work, we are planning to further improve the tsunami simulation code by implementing **GPU computing for several domains to achieve high-resolution**.

# Tsunami Simulation Example
The numerical simulation test defined for an area of 324 x 881 km and 3000 time steps (10min of tsunami propagation) took around 3 minutes of computation.
By modifying the number of time steps in gWave-CPU_mp.f90 (e.g. setting KL parameter to 36000) we can obtain the following inundation for larger area (Peru):

<img width="2360" height="1216" alt="image" src="https://github.com/user-attachments/assets/ba11b2af-8ee8-4372-b67d-1e6fe60fa19b" />

The earthquake slip scenario is defined in the files:
- A_Data/Jimenez.txt based on Jimenez et al (2013).
- A_Data/CaLiBaHu.txt based on Villegas et al (2016).

The hierarchy of the files is:
- gWave-CPU.exe: main program for parallel tsunami simulation.
- 0_files/ Directory where topographic data, bathymetric data and simulation settings is defined.
- 1_deforms/ Directory that has the seafloor vertical deformation divided into partitions (e.g. 48).
- 2_results/ Directory that will store snapshots of tsunami propagation and inundation.

gWave-CPU is a parallel tsunami simulation code based on Imamura (1995).
It was modified for parallel computing by:
- Carlos Davila (carlos.davila.d@uni.pe)
- Julian Palacios (jpalaciose@uni.pe)
- Fernando Garcia (fgarciab@uni.pe)

# References
- Paper submitted to JDR will be referenced here...
- C. Jimenez et al., “Seismic source of 1746 Callao earthquake from tsunami numerical modeling,” J. Disaster Res., Vol.8 No.2, pp. 266-273, 2013. https://doi.org/10.20965/jdr.2013.p0266
- J. C. Villegas-Lanza et al., “Active tectonics of Peru: Heterogeneous interseismic coupling along the Nazca megathrust, rigid motion of the Peruvian Sliver, and Subandean shortening accommodation,” JGR Solid Earth, Vol.121, pp. 7371-7394, 2016. https://doi.org/10.1002/2016JB013080
- Imamura, F., "Review of tsunami simulation with a finite difference method, in Long-Wave Runup Models", edited by H. Yeh, P. Liu, and C. Synolakis, pp. 25–42, 1995. World Scientific Publishing, Hackensack, N. J.
