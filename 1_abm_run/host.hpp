
#ifndef HOST_HPP
#define HOST_HPP

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class Host{
public:
  enum hState{
    hSusceptible =0 , hInfected, hRecovery, hDead //assume lifelong immunity once infected
  };
private:
  hState h_state;
  double tauH;
  //viral dynamics parameters:
  double t_p_host; 
  double t_0_host; 
  double t_r_host; 
  double chi_host; 
  int inf_start_time_host; //infection start time
  int host_id; //add id's to the hosts 
 
  
  
public:
  static int SUSCEPTIBLE_COUNTER;
  static int INFECTIOUS_COUNTER;
  static int RECOVERY_COUNTER;
  static int DEAD_COUNTER; 
  static double ct_lod; 
  static int HOST_ID; //static counter for the hosts 
  
  //viral load model parameters hierarchical parameters for the hosts :
  static double t_p_mean_pop_host; //covariate level time of peak 
  static double t_p_sigma_pop_host; 
  static double t_r_mean_pop_host; 
  static double t_r_sigma_pop_host; 
  static double t_0_mean_pop_host; 
  static double t_0_sigma_pop_host;  
  static double chi_mean_pop_host;
  static double chi_sigma_pop_host; 
  static NumericVector ind_sigamas_host; //individual level variation 
  //add these as
  static NumericVector beta_sigmas_hosts; //sigma to generate individual effect size parameters per viral_kinetic parameters(rows) 
  static NumericVector beta_hosts; //effect size 

  static double t_p_mean_host; //covariate level time of peak 
  static double chi_mean_host; //absolute difference between lod and peak ct
  static double t_r_mean_host; 
  static double t_0_mean_host; 
  
  //add a function to set the viral dynamics parameters to the models
  static void set_viral_static_params_hosts(NumericVector viral_pop_params,
                                      NumericVector ind_sigmas,NumericMatrix betas);
  
  Host(); //create a constructor 
  void set_host_inf_state(Host::hState x_state, int current_time);
  Host::hState get_host_inf_state(); 
  bool check_host_sus();
  bool check_host_inf();
  bool check_host_rec(); 
  bool check_host_dead(); 
  void set_infectious_period_host(double current_time);
  double get_infectious_period_host();
  double get_infection_start_time_host(); 
  void set_host_ct_params();
  //events:
  void infect_host(int current_time);
  void recover_host(double _current_ct);//that is go back to being life long immunity 
  void death_host(); 
  double get_current_ct_host(int current_time); 
  NumericVector get_host_ct_params(); 
  int get_host_id(); 

  //not considering host births and deaths 
};




#endif
