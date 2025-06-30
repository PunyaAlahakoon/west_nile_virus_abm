

 
 ## This folder contains all the codes required to run the ABM with an example. 

## Steps: 
 **_1_abm_structure_:**  
  This folder contains the code for the agent-based model written in Rcpp.  
**Run:**  
  a) **Open** _abm_run.Rproj_ R project.    
  b) **Run** _required_r_package.R_ to install or load the packages required.  
  c) **Run** _run_abm.R_ to run the agent-based model.  
    This script includes and an example based on the parameters used in the model. 


### More details:  
1. **host.hpp:** header declaratons including class and function definitions of the bird agents.
2. **host.cpp:** implementation details of the bird agents.
3. **mosquito.hpp:** header declarations including class and function definitions of the mosquito agents. 
4. **mosquito.cpp:** implementation details of the mosquito agents.
5. **interface.cpp:** implementation of the agent-based model.
    This script calls files 1-4 above to call the characteristics of the bird and mosquito agents.
   Transmission dynamics between mosquitoes and birds, and the mosquitoe life cycle dynamics are captured in this script. The pseudocode can be found in the Supplementary Material.  
   Through this fuction, the required synthetic outputs are called and stored.
6. **run_abm.R:** R script that calls the interface.cpp script. This script runs the agent-based model based on the parameters (as inputs to the interface.cpp) defined in this script.  
   **To run the example,** run this script without making any changes to the model.  
   **To run different simulations,** change the parameters in this script as desired. 
   
