ZZAnalysis: Plotter
===================

To run the plotter, it is enough to run `run_plotter` after having built the code.
Define which plots need to be created and which years you want to plot in the
`run_plotter.cpp`. You can also specify the `Blinded` or `Unblinded` options.

Few important notes:
 * In `Plotter::MakeHistogramsZX` you have to define the `_fs_ROS_SS` and `cb_SS` 
   `std::vectors` needed for the computation of ZX yields. The first vector refers
   to OS/SS ratio and the second one refers to the ratio between the total ZX yield
   and the SS one. These values have to be consistent with the ones we have in `ZXVariables.h` 
   (i.e. update here and there *every time* new ZX yields are computed). If these two
   vectors are properly set, the ZX yields you get from the printout of `Plotter::MakeM4lZX`
   should be the same as the SS+OS combination computed with the scripts in `ZpXEstimation`.
 * To include a new process, define the samples in `run_plotter` and add the corrseponding
   processes in: `Settings.h`, `Plotter::find_current_process` and in `Histogram.cpp`. More
   precisely in the constructors, under `_s_process`, and define the production modes in
   `_s_production_mode`. Also rememeber to update the methods that create legends, otherwise 
   errors will be produced when the script is launched.
 * The script for the yields computation present in this repository is outdated. To compute yields
 see [the reduced trees repository](https://gitlab.cern.ch/HZZ4l/Datacards13TeV/tree/HIG-19-JES).