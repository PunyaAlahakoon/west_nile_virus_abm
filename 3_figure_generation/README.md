
## This folder contains all the R codes that were used to generate the main figures in the manuscript and those in the Supplementary Material. 

**_3_figure_generation_:**  
 To generate the figures, Run:  
   a) **Open** _figure_generation.Rproj_ R project.  
   b) **Run** _required_packages.R_ to load or install the required packages.   
   c) **Run**  _figure_x.R_ (x = 1:7) to generate the figure (number in the same order as in the manuscript) required.  
   d) **Run** _supplement_figures.R_ to generate figures in the Supplementary Material.

     
### More details:  
1. **data folder** contains Bebraska and colorado datasets as well as human cases downloaded from teh CDC website.
2. **functions folder** contains functions to,  
   a) Clean Colorado data _(clean_colorado_data.R)_    
   b) Clean Nebraska data _(clean_nebraska_data.R)_    
   c) ggplot theme to generate figures _(default_theme.R)_  
   d) Approximate the KDE distribution _(kde_approx.R)_  
   e) Calculate the likelihood function to estimate the prevelance _(prev_likelihood_real.R)_  
   f) Sample pools and calculate the pooled ct values _(sample_pools_and_calc_cts.R)_

3. **pre_calculations** folder contains to save time, (but the relevant scripts are still available), pre-calculated values such as kdes and pooled ct values.
4. **synthetic_data** folder contains the synthetic data generated from the ABM (from the example in the run_abm.R file) is in the synthetic data folder.
    
