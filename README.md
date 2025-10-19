# gWave-CPU
CPU-based parallel tsunami simulation code for a single domain.

The numerical simulation test is defined for an area of 324 x 881 km and 3000 time steps (10min of tsunami propagation) took around 3 minutes of computation.
By modifying the number of time steps in 1_7_tun_mp.f90 (e.g. setting KL parameter to 36000) we can obtain the following inundation for Peru:

<img width="1457" height="730" alt="Captura de pantalla 2025-10-19 191811" src="https://github.com/user-attachments/assets/ebcbbd8c-73b0-4150-8303-43da89de73db" />

The earthquake slip scenario is defined in the files:
- A_Data/Jimenez.txt based on Jimenez et al 2013.
- A_Data/CaLiBaHu.txt based on Villegas et al 2016.

The hierarchy of the files is:
- 1_7_tun.exe: main program for parallel tsunami simulation.
- 0_files/ ...
- 1_deforms/ ...
- 2_results/ ...

# How to use
1. Install Fortran compiler with MPI (it can be downloaded from OneAPI)
2. Compile with: mpif90 1_7_tun_mp.f90 1_7_tun_ut.f90 -O3 -o 1_7_tun.exe
3. For the simulation test, run with 48 parallel processes: mpirun -n 48 ./1_7_tun.exe
4. The propagation results will be written in 2_results folder.

Note that the program is only defined for a single domain. For future work, we are planning to further improve the tsunami code by implementing GPU computing and for several domains to achieve high-resolution.

# References
- Paper submitted to JDR will be referenced here...
- Jimenez et al 2013...
- Villegas et al 2016...
- TUNAMI-N2 ...
