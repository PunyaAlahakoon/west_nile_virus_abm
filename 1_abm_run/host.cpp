#include <Rcpp.h>
#include "host.hpp"
//#include <iostream>
using namespace Rcpp;
using namespace sugar;
//using namespace std;

//static variables:
 int Host::SUSCEPTIBLE_COUNTER=0;
 int Host::INFECTIOUS_COUNTER=0;
 int Host::RECOVERY_COUNTER=0; 
 double Host::ct_lod=40; 
 int Host::HOST_ID=0; //initilalise with zero for the host counter 
 
 double Host::t_p_mean_pop_host= log(2.5); 
 double Host::t_p_sigma_pop_host=1; 
 double Host::t_0_sigma_pop_host=0; 
 double Host::t_r_mean_pop_host=log(5);
 double Host::t_r_sigma_pop_host=1; 
 
 NumericVector Host::ind_sigamas_host={.1,.1,.1}; 
 double Host::chi_mean_pop_host=log(12);
 double Host::chi_sigma_pop_host=0.001;
 
 // at the moment, not sampling from the hierarchical model, just assuming that the pop means are the same as the means 
 double Host::t_0_mean_host=0; 
 double Host::t_p_mean_host=t_p_mean_pop_host;
 double Host::t_r_mean_host=t_r_mean_pop_host; 
 double Host::chi_mean_host=chi_mean_pop_host; 
 NumericVector Host::beta_hosts(3,1.1); 
 
 //write a function to set the parameters for the viral dynamics. at the moment, not really setting them 
 
 
//default state of a host/ constructor 
Host::Host(){
  h_state=Host::hState::hSusceptible;
  SUSCEPTIBLE_COUNTER++;
  tauH=0;
  //add the host id:
  host_id=HOST_ID++; 
}

void Host::set_host_inf_state(Host::hState x_state, int current_time){
  h_state=x_state;
  if(x_state==Host::hState::hInfected){
    inf_start_time_host=current_time;  
  }
}

Host::hState Host::get_host_inf_state(){
  return h_state;
}

bool Host::check_host_inf(){
  bool inf_H=false;
  if(h_state==Host::hState::hInfected){
    inf_H=true;
  }
  return inf_H;
}


bool Host::check_host_sus(){
  bool sus_H=false;
  if(h_state==Host::hState::hSusceptible){
    sus_H=true;
  }
  return sus_H;
}

bool Host::check_host_rec(){
  bool inf_R=false;
  if(h_state==Host::hState::hRecovery){
    inf_R=true;
  }
  return inf_R;
}

bool Host::check_host_dead(){
  bool dead_H=false; 
  if(h_state==Host::hState::hDead){
    dead_H=true;
  }
  return dead_H; 
}

void Host::infect_host(int current_time){
  if(h_state==Host::hState::hSusceptible){
    h_state=Host::hState::hInfected;
    inf_start_time_host=current_time;  
    SUSCEPTIBLE_COUNTER--;
    INFECTIOUS_COUNTER++;
  }
}

void Host::recover_host(double _current_ct){
  if((h_state==Host::hState::hInfected) && (_current_ct==ct_lod)){
    h_state=Host::hState::hRecovery;
    INFECTIOUS_COUNTER--;
    RECOVERY_COUNTER++;
  }
}

void Host::death_host(){
  if(h_state==Host::hState::hInfected){
    INFECTIOUS_COUNTER--;
  }else if(h_state==Host::hState::hRecovery){
    RECOVERY_COUNTER--;
  } else if(h_state==Host::hState::hSusceptible){
    SUSCEPTIBLE_COUNTER--; 
      }
  h_state=Host::hState::hDead; //also change the state as well 
}

double Host::get_infectious_period_host(){
  //double tauH=0;
 // if(h_state==Host::hState::hInfected){//randomly create a value from a gamma distrubution 
  // tauH=R::rgamma(2,3);
 // }
  return tauH;
}
  
void Host::set_infectious_period_host(double current_time){
 // double tauH=0;
 // if(h_state==host::hState::hInfected){
    tauH=R::rgamma(2.2,2.5)+current_time; //assume that on average the the infectious perios is 5.5 days 
  //}
}

void Host::set_host_ct_params(){
  NumericVector beta_sig=ind_sigamas_host; 
  //get the parameters for the model:
  t_p_host=abs(R::rlnorm(t_p_mean_host,ind_sigamas_host[0]))+inf_start_time_host;
  chi_host=abs(R::rlnorm(chi_mean_host,ind_sigamas_host[1]));
  t_r_host=abs(R::rlnorm(t_r_mean_host,ind_sigamas_host[2]))+inf_start_time_host;
  /*
  if (chi_host>ct_lod){
    chi_host=exp(chi_mean_pop_host); 
  }
   */
}


double Host::get_current_ct_host(int current_time){
  double ct_current=ct_lod; 
  //check if the host is infected, otherwise the ct =lod; 
  if(h_state==Host::hState::hInfected){

    //minimum ct value to help calculate the slope of the viral load:
    double min_ct= ct_lod-chi_host;

    t_0_host=inf_start_time_host; 
    double omega_p_host=t_p_host-t_0_host; 
    double omega_r_host=t_r_host-t_p_host; 
    
    //find the ct values with time: 
      if(current_time==t_p_host){
      ct_current=min_ct; 
    } else if((current_time>=t_0_host) && (current_time<t_p_host)){
      ct_current=((-chi_host/omega_p_host)*(current_time)) +ct_lod+((chi_host/omega_p_host)*t_0_host);//viral load increasing
    } else if((current_time>t_p_host) && (current_time<=t_r_host)){
      ct_current=((chi_host/omega_r_host)*(current_time))+min_ct-((chi_host/omega_r_host)*t_p_host); 
    }
  }
  
  if(ct_current>ct_lod){
    ct_current=ct_lod; 
  }
   
 return ct_current; 
}

double  Host::get_infection_start_time_host(){
  return inf_start_time_host; 
}

NumericVector Host::get_host_ct_params(){
  NumericVector ct_params(5); 
  ct_params[0]=inf_start_time_host; 
  ct_params[1]=t_0_host; 
  ct_params[2]=t_p_host; 
  ct_params[3]=chi_host; 
  ct_params[4]=t_r_host; 
  
  return ct_params; 
}

int Host::get_host_id(){
  return host_id; 
}
