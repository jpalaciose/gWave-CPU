# gWave-CPU
CPU-based parallel tsunami simulation code for a single domain.

The numerical simulation test is defined for an area of 324 x 881 km and 3000 time steps (10min of tsunami propagation) took around 3 minutes of computation.
By modifying the number of time steps in 1_7_tun_mp.f90 (e.g. setting KL parameter to 36000) we can obtain the following inundation for Peru:

<img width="1457" height="730" alt="Captura de pantalla 2025-10-19 191811" src="https://github.com/user-attachments/assets/ebcbbd8c-73b0-4150-8303-43da89de73db" />

The earthquake slip scenario is defined in the files:
- A_Data/Jimenez.txt based on Jimenez et al 2013.
- A_Data/CaLiBaHu.txt based on Villegas et al 2016.

# References
- Paper submitted to JDR will be referenced here...
- Jimenez et al 2013...
- Villegas et al 2016...
- TUNAMI-N2 ...
