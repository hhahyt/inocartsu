\2D_irregular_waves_bar_break

Included files:
batch.dec - batch run decision file.  If the user wants to run the simulation with no screen output and without running the user interface, change the first number in this file from "0" to "1".  When run in batch mode, the simulation uses the simulation parameters directly from sim_set.dat

bath_loc.m - a Matlab pre-processing script that is used to generate bathymetry input files for COULWAVE from an external source.  See comments inside this file

f_topo.dat - an output file from bath_loc.m; read by COULWAVE to read model bathymetry

nested_grid.grd - an external depth data file that is loaded by bath_loc.m to create the required bathymetry input files for COULWAVE

seqplot_2D.m - a Matlab script file that will plot the model free surface output.  See comments inside this file

sim_set.dat - the primary simulation parameter input file

size_topo.dat - an output file from bath_loc.m; read by COULWAVE to read model bathymetry

spectrum.dat - an input file that provides the discrete amplitude spectrum information required by COULWAVE to run a random wave simulation.  Created with spectrum_2D.m

spectrum_2D.m - a Matlab script file that can be used to generate the spectral inputs required to run a random wave simulation.  This script will generate a 2D TMA spectrum.  See comments inside this file

ts_locations.dat - an input file that stores the spatial locations of to be recorded time series.

x_topo.dat - an output file from bath_loc.m; read by COULWAVE to read model bathymetry

y_topo.dat - an output file from bath_loc.m; read by COULWAVE to read model bathymetry