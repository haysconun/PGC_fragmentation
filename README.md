# PGC_fragmentation

The code in the 'simulation' directory is used to simulate the fragmentation model proposed in 'A fragmentation model for the formation of Drosophila primordial germ cells'.

The function for computing the regular tesselation of the germ plasm gradient is generate_hexagon_tesselation.
Its primary outputs are cell arrays.
cell_output_NC_9_10_11_gp stores the amount of germ plasm associated with each unit (nucleus or cell) over the lineage from NC 9 to NC 11.
cell_output_NC_9_10_11_gp{i,1} stores the NC 9 value for nucleus i.
cell_output_NC_9_10_11_gp{i,2} stores the NC 10 value for each daughter of nucleus i.
cell_output_NC_9_10_11_gp{i,3} store the corresponding values for four offspring at NC 11. 1:2 and 3:4 are sisters, respectively.

cell_output_NC_9_10_11 stores the positions of each unit (nucleus or cell) over the lineage from NC 9 to NC 11.
It is indexed and sorted as cell_output_NC_9_10_11_gp is.

odefun_cell_cycle_oscillotile computes the time derivative of the inhibitory protein concentration and of the cycle phase.
In this function, we specified excess dynamical variables associated with each NC 9 nucleus. 
The dynamics of these excess variables are sequentially activated upon crossings of 2 pi.

solve_embryonic_cell_cycle_oscillator_varying_updated computes the numerical simulation of the growing system of ODEs.
It takes as input precomputed cell_output_NC_9_10_11, cell_output_NC_9_10_11_gp.

sim_examples is an umbrella function that calls generate_hexagon_tesselation and solve_embryonic_cell_cycle_oscillator_varying_updated.
It generates random samples of the entire process from NC 9 to NC 14.

stats_of_cleavage is a function that samples different random displacements of the lattice of nuclei relative to the germ plasm gradient.
It generates statistics of division asymmetries and germ plasm distribution per cell (at NC 11).
It calls generate_hexagon_tesselation but does not explicitely solve the growing system of ODEs.

distinguishable_colors is a helper function that generates a fixed number of distinguishable colors.

plotting_functions/plot_cycle_traces plots dynamics of cycle phases and checks for the number of cycles completed after NC 11.
