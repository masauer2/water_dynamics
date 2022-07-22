# Water Dynamics

## index_waters.c
Identify the H<sub>3</sub>O molecule and all of the H<sub>2</sub>O molecules for the provided simulation frames

## code usage

Run with the following command. Make sure that `id_hydronium.dat` and `water_scan-pos-1.xyz` are included in the same directory.

`./index_waters begin_frame end_frame`

For example, if we want to look at the first 5 frames of the simulation, we would run the following.

`./index_waters 0 5`
