#include <Rcpp.h>
#include "mosquito.hpp"
//#include <iostream>
using namespace Rcpp;
using namespace sugar;
//using namespace std;


//static variables/counters 
int Mosquito::MOZ_ID=0; 
int Mosquito::FEMALE_MOZZIE_COUNTER=0; 
IntegerVector Mosquito::FEAMLE_MOZZIE_IDS; 
int Mosquito::SUSCEPTIBLE_COUNTER=0;
int Mosquito::INFECTIOUS_COUNTER=0;
int Mosquito::DEATHS_COUNTER=0;
//int BIRTH_COUNTER=0; 
double Mosquito::ct_lod=40; 



 double Mosquito::t_p_mean_pop=log(5); // time of peak, log(5)

 double Mosquito::t_p_sigma_pop=1; 
 double Mosquito::t_0_mean_pop=0;
 double Mosquito::t_0_sigma_pop=1; 
 double Mosquito::t_r_mean_pop=log(8+5); //clearance time , log(8+peak). this is not really  being used 
 double Mosquito::t_r_sigma_pop=1; 


 NumericVector Mosquito::ind_sigamas={.1,.2,.3}; 
 double Mosquito::chi_mean_pop=log(15); //absolute difference between lod and peak ct, log(35)
 double Mosquito::chi_sigma_pop=0.001;
 

 //set the parameters at the mean-level 
 //for the moment just taking the absolute values
 double Mosquito::t_0_mean=0; 
 double Mosquito::t_p_mean=t_p_mean_pop;
 double Mosquito::t_r_mean=t_r_mean_pop;
 double Mosquito::chi_mean=chi_mean_pop;

 
 //add betas here:

 NumericVector Mosquito::beta_tarsalis(3,1.1);//t_p, chi, t_r only using this species at the moment 
 NumericVector Mosquito::beta_quinquefasciatus(3,0.025); 
 NumericVector Mosquito::beta_pipiens(3,.025); 
 NumericVector Mosquito::beta_aegypti(3,.025); 
 
 void Mosquito::set_viral_static_params(NumericVector viral_pop_params,
                              NumericVector _ind_sigmas,NumericMatrix _betas){
   t_p_mean_pop=viral_pop_params[0]; // time of peak, log(5)
    t_p_sigma_pop=viral_pop_params[1]; 
    t_0_mean_pop=viral_pop_params[2];
   // t_0_sigma_pop=1; //why did you use this?! 
    t_r_mean_pop=viral_pop_params[3]; //clearance time , log(8+peak)
    t_r_sigma_pop=viral_pop_params[4]; 
    chi_mean_pop=viral_pop_params[5];
    chi_sigma_pop=viral_pop_params[6];
    t_p_mean=viral_pop_params[7];
    t_r_mean=viral_pop_params[8];
    chi_mean=viral_pop_params[9];
    ind_sigamas=_ind_sigmas;
    beta_tarsalis=_betas(0,_); 
    beta_quinquefasciatus=_betas(1,_);
   beta_pipiens=_betas(2,_);
    beta_aegypti=_betas(3,_);
     
 }
 
 
//default state of a mozzie for initial state: no parents are considered for the initial case. 
Mosquito::Mosquito(double prob_events){
  int rand=R::rbinom(1,prob_events);
  //NumericVector rand=Rcpp::runif(1,0,1);
  //Rcpp::Rcout << "Prob events: " << prob_events << std::endl;
  if(rand==1){
    //Rcpp::Rcout <<"Always here" << std::endl;
    ct_bird_change=1;
  }else{
    //Rcpp::Rcout <<"Never here" << std::endl;
    ct_bird_change=2; 
  }
  //bite_digest_p=0; 
  //mozzie_hybernate_period=0; //set the overwitering period to zero 
  //add the age; 
 NumericVector mox_age=Rcpp::rgamma(1,2,3);//mu=6 days, 
 moz_age=floor(mox_age[0]);
  moz_id=MOZ_ID++; //increase the mozzie id every time you create a mozzie object
  moz_parent_id=moz_id;//mozzie id is the same as their id 
  //Rcpp::Rcout << "Creating seed mosquito " << moz_id << std::endl;
  SUSCEPTIBLE_COUNTER++;
  m_state=Mosquito::mState::mSusceptible;
  bird_ct_bitten=ct_lod;
  ct_current=ct_lod;
  overwinter_age=0; 
  NumericVector overwinter_death_agex=Rcpp::rpois(1,114); 
  overwinter_death_age=overwinter_death_agex[0];
  //bernoulli trial for the sex of the mozzie:
 // int bern=R::rbinom(1,0.5); 
  int bern=1; //at the moment, only consider female mosquitoes
  if(bern==1){
  moz_sex=Mosquito::mSex::mFemale;
  FEMALE_MOZZIE_COUNTER++;
  FEAMLE_MOZZIE_IDS.push_back(moz_id);}
  else{
    moz_sex=Mosquito::mSex::mMale; 
  }
  //randomly create the day the mozzie will die: female mozzies will 
  //using a gamma distribution for lifetime 
 // NumericVector death_day=Rcpp::rgamma(1,14,2);//mu=28 days, 
 //NumericVector death_day=Rcpp::rgamma(1,10,3.4);//mu=30.8 days, but skewed distribution
 NumericVector death_day=Rcpp::rgamma(1,5,9);//mu=30.8 days, but skewed distribution
  // death_day=abs(Rcpp::rnorm(1,_avg_lifespan,2)); 
  moz_age_to_die=death_day[0];
  //randomly create a specie out of the four:
  //default species:
  moz_spec=Mosquito::mSpecies::pipiens;
//  moz_spec=static_cast<Mosquito::mSpecies>(rand()%4);
}

 //this is the constructor I use for seeding new infectious mozzies to the population.  
 Mosquito::Mosquito(Mosquito::mState _mstate,int _current_time, double prob_events){//create a female mosqiuo with a given state 
   int rand=R::rbinom(1,prob_events);
   
   //NumericVector rand=Rcpp::runif(1,0,1);
   //Rcpp::Rcout << "Prob events: " << prob_events << std::endl;
   if(rand==1){
     //Rcpp::Rcout <<"Always here" << std::endl;
     ct_bird_change=1;
   }else{
     //Rcpp::Rcout <<"Never here" << std::endl;
     ct_bird_change=2; 
   }
  // bite_digest_p=0; 
//  mozzie_hybernate_period=0; //set the overwitering period to zero 
   bird_ct_bitten=ct_lod;//just initialize with a ct value 
   ct_current=ct_lod;
   overwinter_age=0; 
   NumericVector overwinter_death_agex=Rcpp::rpois(1,114); 
   overwinter_death_age=overwinter_death_agex[0];
   //add the age; 
   moz_age=1; 
   moz_id=MOZ_ID++; 
   moz_parent_id=moz_id;//mozzie id is the same as their id 
   moz_sex=Mosquito::mSex::mFemale;
   FEMALE_MOZZIE_COUNTER++;
   FEAMLE_MOZZIE_IDS.push_back(moz_id);
   m_state=_mstate;
   moz_spec=Mosquito::mSpecies::pipiens;  //default species, if not static cast as above to randomly allocate 
   
   //if the mozzie is infected, set the infected time as well
   if(_mstate==Mosquito::mState::mInfected){
     inf_start_time=_current_time; 
   }
     
 } 
 
 //use this function when mozzies are born. CHANGED:: mozzies goes to the egg state before joining the population
Mosquito::Mosquito(int parent_id, Mosquito::mSpecies parent_sp,Mosquito::mState _mstate,
                   double _vertical_transms_prob,int _current_time, double prob_events ){
  int rand=R::rbinom(1,prob_events);
  
  if(rand==1){
    //Rcpp::Rcout <<"Always here" << std::endl;
    ct_bird_change=1;
  }else{
    //Rcpp::Rcout <<"Never here" << std::endl;
    ct_bird_change=2; 
  }
  bird_ct_bitten=ct_lod;//just initialize with a ct value 
  ct_current=ct_lod;
  overwinter_age=0; 
  NumericVector overwinter_death_agex=Rcpp::rpois(1,114); 
  overwinter_death_age=overwinter_death_agex[0];
  //add the age; 
  moz_age=1; 
  moz_id=MOZ_ID++; //increase the mozzie id every time you create a mozzie object
 if(_mstate==Mosquito::mState::mInfected){
   //run a bernoulli trial 
   int bern=R::rbinom(1,_vertical_transms_prob); // at the moment, there is no vertical transission. howeer, 
   //if decided to include vertical tranmssion, send these as well to the egg state before making them infectious 
   if(bern==1){
     m_state=Mosquito::mState::mInfected; 
     inf_start_time=_current_time; 
     INFECTIOUS_COUNTER++; 
   }else{
     m_state=Mosquito::mState::mEgg;
     
   }
 }
  moz_sex=Mosquito::mSex::mMale;
  //randomly create the day the mozzie will die: male mozzies will 
  //live for 14 days on average 
  NumericVector death_day=Rcpp::rgamma(1,12,1);//mu=12 days
  moz_age_to_die=death_day[0];
  int btr=1; //all the mozzies are female at the moment. 
  if(btr==1){
    moz_sex=Mosquito::mSex::mFemale;
    FEMALE_MOZZIE_COUNTER++;
    FEAMLE_MOZZIE_IDS.push_back(moz_id);
    //randomly create the day the mozzie will die: female mozzies will 
    //live for 21 days on average 
    NumericVector death_day=Rcpp::rgamma(1,5,9);//mu=28 days, 
    moz_age_to_die=death_day[0];
  }
  //mozzies' parent is the same as the new mozzie 
  moz_spec=parent_sp;
  moz_parent_id=parent_id;  //also store the parent's id just in case this is needed 
}
 
void Mosquito::set_egg_period(double birth_time, double rate){ //current_time===birth time 
   NumericVector ee_days=Rcpp::rpois(1,rate); //generate a random day to enter the adult period 
  double e_days=ee_days[0]; 
  double adult_day=birth_time+e_days; //this is the day the 
  egg_p=adult_day; 
 }
 
 double Mosquito::get_egg_time(){
   return egg_p; 
 }
 
 void Mosquito::adult_moz(double current_time, Mosquito::mState _m_state){
   if(m_state==Mosquito::mState::mEgg & egg_p>=current_time){
     m_state=_m_state; 
   }
 }
 
 void Mosquito::set_mozzie_sex_to_female(){
    if(moz_sex==Mosquito::mSex::mMale){
      moz_sex=Mosquito::mSex::mFemale;
      FEMALE_MOZZIE_COUNTER++;
      FEAMLE_MOZZIE_IDS.push_back(moz_id);
    }
 }

void Mosquito::increase_mozzie_age(){
 moz_age++; 
// moz_age +=1; 
}

int Mosquito::get_mozzie_age(){
  return moz_age; 
}

int Mosquito::get_mozzie_death_age(){
  return moz_age_to_die; 
}

Mosquito::mSex Mosquito::get_mozzie_sex(){
  return moz_sex;
}

int Mosquito::get_mozzie_id(){
  return moz_id;
}

Mosquito::mSpecies Mosquito::get_mozzie_species(){
  return moz_spec;
}
  
void Mosquito::set_moz_inf_state(Mosquito::mState y_state){
  m_state=y_state;
}


Mosquito::mState Mosquito::get_moz_inf_state(){
  return m_state;
}

bool Mosquito::check_moz_inf(){
  bool inf_m=false;
  if(m_state==Mosquito::mState::mInfected |m_state==Mosquito::mState::bitten_digesting_infected ){
    inf_m=true;
  } 
  return inf_m;
}

bool Mosquito::check_moz_sus(){
  bool sus_m=false;
  if(m_state==Mosquito::mState::mSusceptible){
    sus_m=true;
  }
  return sus_m;
}

bool Mosquito::check_moz_death(){
  bool death_m=false;
  if(m_state==Mosquito::mState::mDead){
     death_m=true;
  }
  return death_m;
}


void Mosquito::infect_moz(){
  if(m_state==Mosquito::mState::mSusceptible){
     m_state=Mosquito::mState::mInfected;
    SUSCEPTIBLE_COUNTER--;
    INFECTIOUS_COUNTER++; 
    //set the infection start day for the mozzie:
//inf_start_time=_inf_start_time;
  }
}

 void Mosquito::set_inf_start_time(int _inf_start_time){
   if(m_state==Mosquito::mState::mInfected){
   //set the infection start day for the mozzie:
    inf_start_time=_inf_start_time;
    }
 }
 
int Mosquito::get_inf_start_time(){
  return inf_start_time; 
}
 
 //set the mozzie:
 void Mosquito::mozzie_bite_digesting(Mosquito::mState sinfcted, int _current_time){
   if(sinfcted==Mosquito::mState::mInfected){
    m_state= Mosquito::mState::bitten_digesting_infected;
   }if(sinfcted==Mosquito::mState::mSusceptible){
     m_state= Mosquito::mState::bitten_digesting_susceptible; 
   }
 
   //set the digestive period:
  NumericVector digest_p= Rcpp::rpois(1,4);
  //NumericVector xx={5,6,7};
  //NumericVector digest_p= Rcpp::sample(xx,1);
   bite_digest_p=digest_p[0]+_current_time; //add the current time 
   //if(m_state==Mosquito::mState::mInfected){
     //set the infection start day for the mozzie:
   //  inf_start_time=_current_time;s
  // }
 }
 
double Mosquito::get_mozzie_digestive_period(){
  return bite_digest_p; 
}

void Mosquito::recover_moz(){
  if(m_state==Mosquito::mState::mInfected){
    m_state=Mosquito::mState::mSusceptible;
    INFECTIOUS_COUNTER--;
    SUSCEPTIBLE_COUNTER++;
  }
  
}

void Mosquito::death_moz(){
  if(m_state==Mosquito::mState::mSusceptible){
    SUSCEPTIBLE_COUNTER--;
    m_state=Mosquito::mState::mDead; 
    DEATHS_COUNTER++;
  } else if(m_state==Mosquito::mState::mInfected){
    INFECTIOUS_COUNTER--;
    m_state=Mosquito::mState::mDead; 
    DEATHS_COUNTER++;
  }
  if(moz_sex==Mosquito::mSex::mFemale){
    FEMALE_MOZZIE_COUNTER--;
   // int* aa=find(Mosquito::FEAMLE_MOZZIE_IDS.begin(),Mosquito::FEAMLE_MOZZIE_IDS.end(),moz_id);
    //FEAMLE_MOZZIE_IDS.erase(aa); //erase the dead female mozzies 
  }
}

 void Mosquito::set_ct_params(){
     //calculate the parameters:
     NumericVector beta_sig; 
     if(moz_spec==Mosquito::mSpecies::tarsalis){
       beta_sig=beta_tarsalis;
     } else if(moz_spec==Mosquito::mSpecies::quinquefasciatus){
       beta_sig=beta_quinquefasciatus;
     } else if(moz_spec==Mosquito::mSpecies::pipiens){
       beta_sig=beta_pipiens;
     } else if(moz_spec==Mosquito::mSpecies::aegypti){
       beta_sig=beta_aegypti;
     }
     
 
  t_p=abs(R::rlnorm(t_p_mean,ind_sigamas[0]))+inf_start_time;
     chi=abs(R::rlnorm(chi_mean,ind_sigamas[1]));
     t_r=abs(R::rlnorm(t_r_mean,ind_sigamas[2]))+inf_start_time;
     t_0=inf_start_time; 
 }
 
 
 
 double Mosquito::m1_get_current_viral_load(int current_time,double decay_rate){
  double c_t_current=ct_lod; 
 // double v_l_current=1e2; //limit of detection 
  //if the mozzie is not infected, then the ct_value=lod
if(m_state==Mosquito::mState::mInfected){

   //minimum ct value to help calculate the slope of the viral load:
  double min_ct= ct_lod-chi;
  //calculate the viral load curve:
 
  double omega_p=t_p-t_0; 
  double omega_r=t_r-t_p; 
  
  if((current_time>=t_0) && (current_time<=t_p)){
    c_t_current=((-chi/omega_p)*(current_time)) +ct_lod+((chi/omega_p)*t_0);
 //   v_l_cu rrent= peak_vl; //peak viral load 
    
  } else if(current_time<inf_start_time){
    c_t_current=ct_lod;
   // v_l_current=1e2; 
  }
  else{
   if(decay_rate==0){
     c_t_current=min_ct;  
   }else{
    //if decay rate is non-zero, add decay part:
    c_t_current=(decay_rate*current_time)+(min_ct-decay_rate*t_p); 
   }

  }
}
                          
//NumericVector dat={c_t_current,v_l_current};
ct_current=c_t_current;
  return ct_current;
}

 
 NumericVector Mosquito::m22_get_current_viral_load(int current_time,double bird_ct_bitten, double decay_rate){
   double c_t_new=bird_ct_bitten; 
   double c_t_current; 
   if(decay_rate==0){
     c_t_current=c_t_new; 
   }else{
  //calculate the ct value current 
  c_t_current=(decay_rate*current_time)+c_t_new-(decay_rate*t_0); 
     c_t_current=std::min(c_t_current,ct_lod); 
   }
   
  // c_t_current=35;
  // return c_t_current; 
   double ct_method= 2; 
   return {c_t_current,ct_method};
 //  Rcpp::Rcout << "c_t_current: " << c_t_current << std::endl;
 }
 
 NumericVector Mosquito::m33_get_current_viral_load(double host_current_ct, double prob_events, int current_time,double decay_rate){
   //generate a random number:
   double c_t_current; 
   double ct_method; 

     //Rcpp::Rcout << "Option 2, inherited Ct" <<  std::endl;
     ct_method=3;

     //minimum ct value to help calculate the slope of the viral load:
     double min_ct= host_current_ct-chi;//peak viral load 
     //calculate the viral load curve:
     double omega_p=t_p-t_0; 
     double omega_r=t_r-t_p; 
     
     if((current_time>=t_0) && (current_time<=t_p)){
       c_t_current=((-chi/omega_p)*(current_time)) +host_current_ct+((chi/omega_p)*t_0);
     } else if(current_time<inf_start_time){
       c_t_current=ct_lod;
       // v_l_current=1e2; 
     } else if(current_time==inf_start_time){
       c_t_current=host_current_ct; 
     }
     else{
       double c_t_c=(decay_rate*current_time)+(min_ct-(decay_rate*t_p)); 
       c_t_current=std::min(c_t_c,ct_lod);
       //v_l_current=((peak_vl-1e2)/omega_p)*(current_time)+1e2-((peak_vl-1e2)/omega_p)*(t_0); 
       
     }
 
   ct_current = c_t_current;
   return {ct_current,ct_method};
 }
 
 NumericVector Mosquito::get_mozzie_ct_params(){
   NumericVector ct_params(5);
   ct_params[0]=inf_start_time;
   ct_params[1]=t_0;
   ct_params[2]=t_p;
   ct_params[3]=chi;
   ct_params[4]=t_r; 
   
   return ct_params; 
 }

void Mosquito::increase_mozzie_bite_counter(){
  mozzie_bite_counter++; 
}

int Mosquito::get_mozzie_bite_counter(){
  return mozzie_bite_counter; 
}

 
void Mosquito::set_bitten_bird_ct_to_mozzie(double current_bird_ct,double viral_percent, double prob_events){
  double lg_vl_per;
  if(viral_percent==0){
    lg_vl_per=0;
  }else{
    lg_vl_per=log10(viral_percent);
  }
  double vl_frac = ((current_bird_ct - 36.9)/-2.7) + lg_vl_per;
  double new_ct = 36.9 - 2.7 * vl_frac;
  bird_ct_bitten=new_ct;
  
} 
 
 double Mosquito::get_bitten_bird_ct_to_mozzie(){
   return bird_ct_bitten; 
 }
 
 int Mosquito::get_ct_model_change(){
   //ct_bird_change=2;
   return ct_bird_change; 
  
 }
 
 void Mosquito::mozzie_overwinter(){
   if(m_state==Mosquito::mState::mSusceptible){
     m_state=Mosquito::mState::mOverwintering_susceptible; 
   }else if(m_state==Mosquito::mState::mInfected){
     m_state=Mosquito::mState::mOverwintering_infected; 
   }
 }
 
 
 void Mosquito::set_overwinter_age(int _overwinter_age){
   overwinter_age=_overwinter_age; 
 }
 
 void Mosquito::increase_overwinter_age(){
   overwinter_age++; 
 }
 
 int Mosquito::get_overwinter_age(){
  return overwinter_age; 
 }
 
 
 void Mosquito::mozzie_become_active(){
   if(m_state==Mosquito::mState::mOverwintering_susceptible){
     m_state=Mosquito::mState::mSusceptible;
   }else if(m_state==Mosquito::mState::mOverwintering_infected){
     m_state=Mosquito::mState::mInfected; 
   }
 }
 
 
 int Mosquito::get_overwinter_death_age(){
   return overwinter_death_age; 
 }
 