# Asexual-Gaussian-Selection

The C++ program Gaussian.cpp enables one to evaluate the long-term steady-state distribution of mean phenotypes for an asexual population under Gaussian selection and was used to generate all of the primary results in the paper “The Divergence of Mean Phenotypes Under Persistent Gaussian Selection” by Lynch and Menor. A verbal description of the model can be found at the top of the code. Model parameters for individual runs can be entered at the #define lines at the top of the program. The program is written to simultaneously make parallel runs at 21 different population sizes. 

All of the component parts should be sitting in a single folder.

If one wishes to simultaneously make runs with various combinations of parameters, the program launch.sh can be used by entering ./launch.sh at the command line. The parameters that can be run match those in the #define terms at the top of the program, and desired values need to be entered into the “Parameters” file; when multiple values for a parameter are not to be used, the value entered in the #define line are used as defaults.

Gaussian.sh will need to be edited for one’s particular system.

Modules.sh loads the GSL libraries necessary for running the program.

The output will appear in a folder titled Runs, which gives the final output for a run (dataout), as well as Slurm files which record cumulative estimates over time (enabling the user to determine if runs have been sufficiently long to yield stable results)

Concatenate.sh can be used to produce a single file containing all of the dataout files within a run, which can then be entered into a spreadsheet, etc. 

