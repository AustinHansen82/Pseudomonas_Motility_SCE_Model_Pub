# Pseudomonas_Motility_SCE_Model_Pub
Code for "Modeling Study of the Role of Wrap Mode and Reversals in Pseudomonas aeruginosa Motion"
The main code for modeling the migration of Pseudomonas aeruginosa is written in C++. The code is set up to simulate the migration of 100 bacterium towards a simple point source of chemoattractant. The concentration of the chemoattractant increases linearly in the radial direction toward the source. 

To run the code and reproduce Figure 4ABC of the paper
  1) Download the repository folder
  2) In the folder SCE_Code, compile the code through the command "g++ *.cpp -std=c++11"
  3) Run each simulation by the command "./a.out -slurm #" where # is the simulation number (1 and 2)
  4) Once the two simulations have finished, run the matlab script titled "Bacteria_PostProcess_Combined_Figure4ABC" - this will generate a folder and stored the figure in the folder 

In order for the code to successfully run you must have the following folders in with the code
  1) "animation" - this should contain a list of folders ("machine1, machine2 , ... , ) where each is a folder for simulation (run1, run2 , ... ,) - Each of these folders stores the vtk files to visualize the simulations
  2) "dataStats" - this should contain a list of folders ("machine1, machine2 , ... , ) where each is a folder for simulation (run1, run2 , ... ,) - Each of these folders stores the data files to analyze the simulations
  3) "resources" - this should contain .cfg files (test1.cfg, test2.cfg, ... , ) where each is a parameter file used to run simulation (run1, run2 , ... ,)

Notes:   
  1) There is a folder titled - "Simulation_Resource_Files" which contains all of the parameter files run for the simulations in the paper "Modeling Study of the Role of Wrap Mode and Reversals in Pseudomonas aeruginosa Motion"
  2) There is a folder titled - "Chemoattractant_Concentrations" which contains all of the chemoattractant concentrations used in the paper "Modeling Study of the Role of Wrap Mode and Reversals in Pseudomonas aeruginosa Motion". When using these, the folder location will need to be changed to match your environment. This can be found in the file "TissueGrid.cpp" in the function "Create_Experimental_Gradient"
  3) Some of the simulations in the paper "Modeling Study of the Role of Wrap Mode and Reversals in Pseudomonas aeruginosa Motion" require minor changes to the code based based on the simulation designs. 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14768453.svg)](https://doi.org/10.5281/zenodo.14768453)
