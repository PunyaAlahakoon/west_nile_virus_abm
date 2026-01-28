
#ifndef MOSQUITO_HPP
#define MOSQUITO_HPP

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class Mosquito{
public:
  enum mState{
    mSusceptible =0 , mInfected, mDead, mEgg, mOverwintering_infected, mOverwintering_susceptible, bitten_digesting_infected,bitten_digesting_susceptible
  };
  enum mSex{
    mMale=0, mFemale=1
  };
  enum mSpecies{
    tarsalis=0, quinquefasciatus, pipiens, aegypti 
  };
  
private:
  mState m_state;
  mSex moz_sex;
  mSpecies moz_spec;
  int moz_id; //each mozzie has an id
  double egg_p; //period from egg to adult 
  double bite_digest_p; //period of digesting and laying of eggs 
 // Mosquito* parent;
  int moz_parent_id; 
  int moz_age; //mozzie age 
  int moz_age_to_die; //randomly allocate an age to die as well when creating a mozzie 
  NumericVector viral_dynamics;
  //viral dynamics individual level parameters:
  double t_p; 
 // double omega_p; 
  double t_0; 
  //double omega_r; 
  double t_r; 
  double chi; 
  int inf_start_time; //infection start time
  int mozzie_bite_counter; 
  double bird_ct_bitten; 
   double ct_current; 
   int ct_bird_change; 
   //double mozzie_hybernate_period; 
 // NumericMatrix design_matrix; 
 //initiliase an overwintering age for the mozzies:
 int overwinter_age; 
 int overwinter_death_age; 

public:
  static int MOZ_ID; 
  static int FEMALE_MOZZIE_COUNTER; 
  static IntegerVector FEAMLE_MOZZIE_IDS; 
  static int SUSCEPTIBLE_COUNTER;
  static int INFECTIOUS_COUNTER;
  static int DEATHS_COUNTER;

  //static int BIRTH_COUNTER;
  
  //viral load model parameters hierarchical parameters:
  static double t_p_mean_pop; //covariate level time of peak 
  static double t_p_sigma_pop; 
  static double t_r_mean_pop; 
  static double t_r_sigma_pop; 
  static double t_0_mean_pop; 
  static double t_0_sigma_pop;  
  static double chi_mean_pop;
  static double chi_sigma_pop; 
  static NumericVector ind_sigamas; //individual level variation 
  //add these as
  static NumericVector beta_sigmas_tarsalis; //sigma to generate individual effect size parameters per viral_kinetic parameters(rows) 
  static NumericVector beta_tarsalis; //effect size 
  static NumericVector beta_sigmas_quinquefasciatus;
  static NumericVector beta_quinquefasciatus; 
  static NumericVector beta_sigmas_pipiens;
  static NumericVector beta_pipiens; 
  static NumericVector beta_sigmas_aegypti;
  static NumericVector beta_aegypti; 

  static double t_p_mean; //covariate level time of peak 
  static double chi_mean; //absolute difference between lod and peak ct
  static double t_r_mean; 
  static double t_0_mean; 

    
  static void set_viral_static_params(NumericVector viral_pop_params,
                                      NumericVector ind_sigmas,NumericMatrix betas);



  
  static double ct_lod; 
  //maybe try having a static ct current as well 

  
  
  Mosquito(double prob_events); //create a constructor for initial mozzies who do not have parents from the current pop. 
  Mosquito(Mosquito::mState _mstate,int _current_time, double prob_events); //create a female mosqiuo with a given state 
  Mosquito(int parent_id, Mosquito::mSpecies parent_sp,Mosquito::mState _mstate,
           double _vertical_transms_prob, int _current_time, double prob_events);
  void increase_mozzie_age(); //increase the age by one day; 
  int get_mozzie_age();
  int get_mozzie_death_age(); 
  
  void set_inf_start_time(int _inf_start_time); 
  void mozzie_bite_digesting(Mosquito::mState should_the_mozzie_be_infcted, int _inf_start_time); 
  
  void set_overwinter_age(int _overwinter_age); 
  void increase_overwinter_age(); 
  int get_overwinter_age(); 
  int get_ct_model_change(); 
  int get_female_mozzie_number();
  void set_mozzie_sex();
  Mosquito::mSex get_mozzie_sex();
  int get_mozzie_parent_id(int n_female_mozzies_at_t); //randomly get an id of a female mozzie that can give birth 
  Mosquito::mSpecies get_mozzie_species();
  void set_mozzie_species();
  int get_mozzie_id();
  int get_random_parent_id(); //create a random female mozzie iD to become a parent 
  
  void set_moz_inf_state(Mosquito::mState y_state);
  Mosquito::mState get_moz_inf_state();
  bool check_moz_sus();
  bool check_moz_inf(); 
  bool check_moz_death();
  double get_infectious_period_moz();
  void set_infectious_period_moz(double current_time);
  //void viral_load_curve_moz(IntegerVector ages,NumericVector kinetic_pars,int times);
  void viral_load_curve_moz(int st_time,int times);
  double current_viral_load(int current_time,double host_current_ct, int method); 
  NumericVector m22_get_current_viral_load(int current_time,double bird_ct_bitten,double decay_rate);
  NumericVector m33_get_current_viral_load(double host_current_ct, double prob_events, int current_time, double decay_rate); 
  
  //folllowing there functions are actually not being used at the moment. the above one is being used 
  double m1_get_current_viral_load(int current_time,double decay_rate);//ct model without considering the ct 
  //values of the hosts 
  NumericVector m2_get_current_viral_load(double host_current_ct, double prob_events,int current_time);

  //ct model starts with the bird ct 
 NumericVector get_mozzie_ct_params(); 
  int get_inf_start_time(); 
  void set_mozzie_sex_to_female(); 
  void set_bitten_bird_ct_to_mozzie(double current_bird_ct,double viral_percent, double prob_events); 
  double get_bitten_bird_ct_to_mozzie(); 
 
 int get_overwinter_death_age(); 
  //set the viral dynamics params when infected:
  void set_ct_params();
  //events:
  void infect_moz();
  void recover_moz();//go back to being susceptible 
  void death_moz();
  void adult_moz(Mosquito::mState y_state); //the egg become an adult 
  void increase_mozzie_bite_counter(); 
  int get_mozzie_bite_counter(); 
  double get_egg_time(); //time period for the egges to become adults 
  void set_egg_period(double current_time, double rate); //get the period for the eggs to become adults 
 //int birth_moz();//birth of one mozzie 
  void adult_moz(double current_time,Mosquito::mState _m_state); 
  double get_mozzie_digestive_period(); 
  void mozzie_overwinter();
  //double get_overwintering_period(); 
  void mozzie_become_active(); //end of the overwintering period 
  //make the mozzie susceptible from egg state, also kepe the option to change the 
  //mozzie state in case you need to consider vertical transmission 
  
};



#endif
