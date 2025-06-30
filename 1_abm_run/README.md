

 
 ## This folder contains all the code required to run the ABM with an example. 


 **_1_abm_run_:**  
  This folder contains the code for the agent-based model written in Rcpp.  
**Run:**  
  a) **Open** _abm_run.Rproj_ R project.    
  b) **Run** _required_r_package.R_ to install or load the packages required.  
  c) **Run** _run_abm.R_ to run the agent-based model.  
    This script includes an example based on the parameters used in the model. 


### More details:  
1. **host.hpp:** Header declarations including class and function definitions of the bird agents.
2. **host.cpp:** Implementation details of the bird agents.
3. **mosquito.hpp:** Header declarations including class and function definitions of the mosquito agents. 
4. **mosquito.cpp:** Implementation details of the mosquito agents.
5. **interface.cpp:** Implementation of the agent-based model.
    This script calls files 1-4 above to call the characteristics of the bird and mosquito agents.
   Transmission dynamics between mosquitoes and birds, and the mosquitoe life cycle dynamics are captured in this script. The pseudocode can be found in the Supplementary Material.  
   Through this function, the required synthetic outputs are called and stored.
6. **run_abm.R:** R script that calls the interface.cpp script. This script runs the agent-based model based on the parameters (as inputs to the interface.cpp) defined in this script.  
   **To run the example,** run this script without making any changes to the model.  
   **To run different simulations,** change the parameters in this script as desired. 
   

#### Other folders:  
1. **functions:** This folder contains the sinusoidal functions used for overwintering and other parameters in the ABM.  
2. **mozzie_outputs:** This folder will store the synthetic data related to the mosquitoes after running the ABM:   
   a) **mozzie_object_data:** Stores the mosquito object characteristics at each time-step.  
   b) **params:** Stores the parameters of some of the mosquito objects' Ct model.  
   c) **states:** Store the prevalence, incidence, total population, number of susceptible, deaths, etc, of the mosquito population at each time-step.
3. **host_outputs:** This folder will store the synthetic data related to the birds after running the ABM:   
   a) **host_object_data:** Stores the bird object characteristics at each time-step.  
   b) **params:** Stores the parameters of some of the bird objects' Ct model.  
   c) **states:** Store the prevalence, incidence, total population, number of susceptible, deaths, etc, of the bird population at each time-step.  
