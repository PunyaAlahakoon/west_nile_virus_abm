# Tracking West Nile Virus dynamics using viral loads from trapped mosquitoes

### This folder contains all the codes and synthetic data generated, and observed data used in this manuscript by P. Alahakoon, I. Marchinton, J. Fauver, and J. Hay.  

#### Open the folders in their numbered order. Each folder has a README file that contains more information about the content of the folder. 

### <ins>Folder structure</ins>

1. **_1_abm_run_:**  
  This folder contains the code for the agent-based model written in Rcpp.  
**Run:**  
  a) **Open** _abm_run.Rproj_ R project.    
  b) **Run** _required_r_package.R_ to install or load the packages required.  
  c) **Run** _run_abm.R_ to run the agent-based model.  
    This script includes an example based on the parameters used in the model. 

   
2. **_2_kde_calculations_:**  
  **Run:**  
   a) **Open** _2_kde_calculations.Rproj_ R project.  
   b) **Run** _kde_pre_calculations.R_ to calculate and store KDE approximations based on different pool sizes and the number of positive samples in a pool. 

   
3. **_3_figure_generation_:**  
  **Run:**  
   a) **Open** _figure_generation.Rproj_ R project.  
   b) **Run** _required_packages.R_ to load or install the required packages.   
   c) **Run**  _figure_x.R_ (x = 1:7) to generate the figure (number in the same order as in the manuscript) required.  
   d) **Run** _supplement_figures.R_ to generate figures in the Supplementary Material.

