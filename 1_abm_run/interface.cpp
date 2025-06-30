

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "mosquito.hpp"
#include "host.hpp"
#include <math.h>
//#include "pooled_ct_calc.cpp"

using namespace Rcpp;
using namespace std; 
using namespace sugar;


// [[Rcpp::export]]
void run_mozzie_model(int time, int nhosts, int nMozzies, NumericVector daily_pr_bite,
                      NumericVector mozzie_birth_rate_t,
                      double probhostInfection, NumericVector probMozzInfection_t, 
                      int iniInfectedMozzies,
                      double host_birth_rate,double host_death_rate,
                      double vertical_transms_prob,
                      NumericVector daily_pr_death, 
                      NumericVector seed_on,
                      bool seed_mozzies,
                      int sim_number, double egg_to_adult_days, 
                      int mozzie_viral_dynamics_method,double prob_ct_model_change,
                      double viral_percent_from_birds,NumericVector set_the_overwintering_period,
                      NumericVector mozzie_owerwintering_prob,NumericVector active_period,
                      NumericVector prob_mozzies_become_active,
                      bool VERBOSE){
  
  
  //List outputs=List::create(0);
  if(VERBOSE){
    Rcpp::Rcout << "========== Initial conditions ==========" << std::endl;
    Rcpp::Rcout << "Starting simulation with " << nhosts << " hosts and " << nMozzies << " mosquitoes" << std::endl;
    Rcpp::Rcout << "Initial conditions: " << iniInfectedMozzies << " infected mosquitoes" << std::endl;
    Rcpp::Rcout << "========================================" << std::endl<< std::endl;
  }
  
  
  double dt=1;   //timestep 
  //double alpha=0.9;//amplitude of seasonal variation 
  double t_rate=0; //time dependent birth rate with seasonality 
  double prob_death_hosts=1-exp(-host_death_rate*dt); 
  
  
  if(VERBOSE){
    Rcpp::Rcout << "========== Model parameters ==========" << std::endl;
    Rcpp::Rcout << "Duration: " << time << " time steps" << std::endl;
    //Rcpp::Rcout << "Mosquito daily bite probability: " << prBite << std::endl;
    //  Rcpp::Rcout << "Probability of host infection given bite: " << probhostInfection << std::endl;
    //    Rcpp::Rcout << "Probability of mosquito infection given bite: " << probMozzInfection << std::endl;
    Rcpp::Rcout << "======================================" << std::endl<< std::endl;
    
  }
  
  //create ofstream for file outputs 
  std::ofstream host_states;
  stringstream host_ss;
  host_ss << "host_outputs/states/host_states_"  << sim_number << "_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string file_hosts= host_ss.str();
  host_states.open (file_hosts);
  host_states<<"Susceptible,Infectious,Recovered,Deaths,Births,Total_pop,new_inf_hosts, \n"; 

  std::ofstream mozzie_states; 
  stringstream mozzie_ss;
  mozzie_ss << "mozzie_outputs/states/mozzie_states_"  << sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string file_mozzies=mozzie_ss.str();
  mozzie_states.open(file_mozzies);
  mozzie_states<<"Susceptible, Infectious,Deaths, Births, new_adults,Total_pop,new_inf_mozzies,\n" ; 
  
  std::ofstream mozzie_ct_params; 
  stringstream mozzie_pars; 
  mozzie_pars << "mozzie_outputs/params/mozzie_ct_params_"  << sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string file_mozzie_pars=mozzie_pars.str(); 
  mozzie_ct_params.open(file_mozzie_pars);
  mozzie_ct_params<<"mozzie_id,current_time,inf_start_time,t_0,t_p,chi,t_r,\n"; 
  
  std::ofstream host_ct_params; 
  stringstream host_pars;
  host_pars << "host_outputs/params/host_ct_params_" << sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string file_host_pars=host_pars.str();
  host_ct_params.open(file_host_pars); 
  host_ct_params<<"host_id,current_time,inf_start_time,t_0,t_p,chi,t_r,\n";
  

  
  
  //IntegerVector ages=seq_len(time);
  IntegerMatrix timeStatehosts(time+1,5);  //initialise a matrix to store states of hosts
  //suscetible, infectious, recovered, dead, births 
  IntegerMatrix timeStateMozzies(time+1,5); //for mozzies== susceptible, infectious, dead, births, 
  NumericVector nInfIniM; 
  
  vector<Host*> current_hosts;
  vector<Mosquito*> current_mosquitoes;
  vector<Mosquito*> dead_mosquitoes;
  vector<Mosquito*> egg_mozzies; 
  vector<Mosquito*> ovewintering_mozzies;
  vector<Mosquito*> digesting_mozzies; 
  
  NumericVector infectiousPeriods_mozzie; 
  int mozzie_births; 
  //NumericVector pro_female_mozzies; //number of female mozzies over time
  NumericVector trapped_mozzies;//number of mozzies that get trapped
  NumericVector trapped_mozzie_cts; //ct values of the traped mozzies 
  
  //captured mozzies with their characteristics:
  NumericVector median_ct_moz_tarsalis; 
  NumericVector median_ct_moz_quinquefasciatus;
  NumericVector median_ct_moz_pipiens;
  NumericVector median_ct_moz_aegypti;
  
  NumericVector total_mozzies;
  NumericVector new_infs; 
  
  //would be good to have these in a dataframe 
  NumericVector time_since_infection; 
  NumericVector time_of_ct_collection; 
  
  
  NumericVector death_ages; 
  NumericVector mozzie_ages; 
  nInfIniM=Rcpp::sample(nMozzies,iniInfectedMozzies);//replace=false 
  //set the initial inectious mozzies
  
  Mosquito* new_mozzie;
  for(int j=0;j<nMozzies;j++){
    new_mozzie =new Mosquito(prob_ct_model_change); //create a mozzie without parents. all these are females.
    current_mosquitoes.push_back(new_mozzie);
    if((any(j==nInfIniM)).is_true()){
      //change the state to infectious
      current_mosquitoes[j]->infect_moz();//set the mozzie as infectious and set the infection start day to zero
      current_mosquitoes[j]->set_inf_start_time(0);
      //set the infected mozzie's gender to female in case they were male
    }
  }
  
  // all the hosts are initially susceptible 
  Host* new_host;
  for(int k=0;k<nhosts;k++){
    new_host=new Host();
    current_hosts.push_back(new_host);
    // inihostState[k].set_host_inf_state(host::hState::hSusceptible);
  }
  
  //add the initial numbers to the matrices  
  IntegerVector iniHostN={nhosts,0,0,0,0};
  IntegerVector iniMozN={nMozzies-iniInfectedMozzies,iniInfectedMozzies,0,0,0};
  timeStatehosts(0,_)=iniHostN;
  timeStateMozzies(0,_)=iniMozN;
  //pro_female_mozzies.push_back((float)iniInfectedMozzies/nMozzies);
  
  // add the initial states 
  //  mozzie_states<<"Susceptible, Infectious,Deaths, Births, new_adults, Total_pop,\n" ; 
  mozzie_states<<nMozzies-iniInfectedMozzies << "," << iniInfectedMozzies << "," << 0 << "," << 0 << "," << 0<< ","<< nMozzies << "," << 0 << "\n"; 
  host_states<<nhosts << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << nhosts << "," << 0 << "," <<"\n"; 
  
 // int female_moz_pop=iniInfectedMozzies; 
  total_mozzies.push_back(nMozzies);
  
  //open a file to store mozzie objects and store the data. ditto for hosts 
  std::ofstream mozzie_objects_i; 
  stringstream ss;
  ss << "mozzie_outputs/mozzie_object_data/mozzie_objects_" << 0 << "_"<< sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string filename= ss.str();
  mozzie_objects_i.open(filename); 
  //mozzie_objects_i.open("mozzie_outputs/mozzie_object_data/mozzie_states.csv"); 
  //add column names; 
  mozzie_objects_i << "id,inf_state,current_time,infected_time,ct_value,ct_type_method, number_of_bites,species,sex,current_age,death_day,sim_number,\n"; 
  //get the values needed to store the dataset; 
  for(int m=0;m<current_mosquitoes.size();m++){
    //get the ct value of the bird who bit at that time of bite:
    double ct_bird_moz=current_mosquitoes[m]->get_bitten_bird_ct_to_mozzie(); 
    double mozzie_ct=40; 
    int ct_type=0; 
    //get the current ct value of the mozzie based on the method:
    if(mozzie_viral_dynamics_method==1){
      mozzie_ct=current_mosquitoes[m]->m1_get_current_viral_load(0);
      ct_type=1;
    }else{
      int ct_ch =current_mosquitoes[m]->get_ct_model_change(); 
      //Rcpp::Rcout << "ct_bird_return: " << ct_ch << std::endl;
      NumericVector mozzie_ct_dat; 
      if(ct_ch==1){
        mozzie_ct_dat=current_mosquitoes[m]->m22_get_current_viral_load(ct_bird_moz);
        // mozzie_ct_dat={35,2};
      }else{
        mozzie_ct_dat=current_mosquitoes[m]->m33_get_current_viral_load(ct_bird_moz,prob_ct_model_change,0); 
        // mozzie_ct_dat={35,3};
      }
     // NumericVector mozzie_ct_dat=current_mosquitoes[m]->m2_get_current_viral_load(ct_bird_moz,prob_ct_model_change,0);
      mozzie_ct=mozzie_ct_dat[0];
      ct_type=mozzie_ct_dat[1];
    }
    mozzie_objects_i << current_mosquitoes[m]->get_mozzie_id() << "," <<
      current_mosquitoes[m]->get_moz_inf_state() << "," <<
        0 << "," <<   current_mosquitoes[m]->get_inf_start_time() << "," <<
          mozzie_ct << "," << ct_type <<  "," <<current_mosquitoes[m]->get_mozzie_bite_counter() << "," <<
              current_mosquitoes[m]->get_mozzie_species() << ","  <<
                current_mosquitoes[m]->get_mozzie_sex() << "," <<
                  current_mosquitoes[m]->get_mozzie_age() << "," <<
                  current_mosquitoes[m]->get_mozzie_death_age() << "," << 
                  sim_number << "\n"; 
  }
  mozzie_objects_i.close();
  
  
  std::ofstream host_objects_i; 
  stringstream hh; 
  hh << "host_outputs/host_object_data/host_objects_"  << 0 << "_"<< sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
  string file=hh.str(); 
  host_objects_i.open(file); 
  host_objects_i << "id,inf_state,current_time,infected_time,ct_value,sim_number,\n";
  for(int m=0; m<current_hosts.size();m++){
    host_objects_i<<current_hosts[m]->get_host_id()  << "," << 
      current_hosts[m]->get_host_inf_state() << "," << 0 << "," << 
        current_hosts[m]->get_infection_start_time_host() << "," << 
          current_hosts[m]->get_current_ct_host(0) << "," << 
            sim_number  << "\n"; 
  }
  host_objects_i.close(); 
  
  
  //increment the event with times 
  int i=1; 
  while((i <time+1)){//make sure that there's a non-zero female pop all the time 
    //at each time point, open a file for the mozzie objects:
    std::ofstream mozzie_objects_i; 
    stringstream ss;
    ss << "mozzie_outputs/mozzie_object_data/mozzie_objects_" << i << "_"<< sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
    string filename= ss.str();
    mozzie_objects_i.open(filename); 
    //mozzie_objects_i.open("mozzie_outputs/mozzie_object_data/mozzie_states.csv"); 
    //add column names; 
    mozzie_objects_i << "id,inf_state,current_time,infected_time,ct_value,ct_type_method,number_of_bites,species,sex,current_age,death_day,sim_number,\n" ; 
    //ditto for hosts:
    std::ofstream host_objects_i; 
    stringstream hh; 
    hh << "host_outputs/host_object_data/host_objects_"  << i << "_"<< sim_number <<"_" << prob_ct_model_change << "_" << viral_percent_from_birds << ".csv";
    string file=hh.str(); 
    host_objects_i.open(file); 
    host_objects_i << "id,inf_state,current_time,infected_time,ct_value,sim_number,\n";
    
    int new_inf_mozzies=0; 
    int new_inf_hosts=0; 
    //    Rcpp::Rcout << "Time step: " << i << std::endl;
    //events to consider:
    //mozzie bites host 
   // int nDeaths=0; 
    for(int j=0;j<current_mosquitoes.size();j++){
      NumericVector rands=runif(3,0,1);
    //  if(current_mosquitoes.size()>0){
      if(current_mosquitoes[j]->get_mozzie_sex()==Mosquito::mSex::mFemale){//if the mozzies doesn't die, and the mozzie is male, increase the age
        //if the mozzie is a female, it may bite a host 
        //calculate the biting rate:
        double prBite=daily_pr_bite[i]; 
        if(rands[0]<=prBite) { //the mozzie has bitten the host
          current_mosquitoes[j]->increase_mozzie_bite_counter(); 
          //generate a random host:
          IntegerVector rhost_v=Rcpp::sample(nhosts,1);
          int rhost=rhost_v[0]-1; //make sure to check if the host was included somewhere in the same time step
          //now check if the mozzie is infected or not and the host is infected or not:
          if((current_hosts[rhost]->get_host_inf_state()==Host::hState::hInfected) && 
             (current_mosquitoes[j]->get_moz_inf_state()==Mosquito::mState::mSusceptible)){//if the host is infected and the mozzie is susceptible
            if(rands[1]<=probMozzInfection_t[i]){//infect the mozzie 
              // if(rands[2]<=0.5){//infect the mozzie 
              //if the mozzie is infected, send it to the mozzie digestive period: as 
              current_mosquitoes[j]->infect_moz();
              //set the infection start time:
              current_mosquitoes[j]->set_inf_start_time(i);
            //  current_mosquitoes[j]->infect_moz(i);
              new_inf_mozzies++; 
              //store the bird ct value as a characteristic to the mozzie 
              double bird_ct=current_hosts[rhost]->get_current_ct_host(i); 
              current_mosquitoes[j]->set_bitten_bird_ct_to_mozzie(bird_ct,viral_percent_from_birds,prob_ct_model_change); 
              current_mosquitoes[j]->set_ct_params(); 
            } 
          } else if((current_hosts[rhost]->get_host_inf_state()==Host::hState::hSusceptible) && 
            (current_mosquitoes[j]->get_moz_inf_state()==Mosquito::mState::mInfected)){ //host gets infected
            if(rands[2]<=probhostInfection){
              current_hosts[rhost]->infect_host(i);
              current_hosts[rhost]->set_host_ct_params(); 
              new_inf_hosts++; 
              //attach an infectious period to the host 
              //current_hosts[rhost]->set_infectious_period_host(i);// no need to do this anymore 
            } 
          } 
        //at the end of the succefful bite, remove the mozzie to from the current mozzies and add it to the digesting mozzie vector:
        //set the mozzie as to digesting
        current_mosquitoes[j]->mozzie_bite_digesting(current_mosquitoes[j]->get_moz_inf_state(),i);
        digesting_mozzies.push_back(current_mosquitoes[j]);
        current_mosquitoes.erase(current_mosquitoes.begin()+j); 
        }
      }
    }
    
    //check if the digesting mozzies can go back to the current population:
    for(int a=0; a<digesting_mozzies.size();a++){
      if(digesting_mozzies[a]->get_mozzie_digestive_period()<=i){
        if(digesting_mozzies[a]->get_moz_inf_state()==Mosquito::mState::bitten_digesting_infected){
          digesting_mozzies[a]->set_moz_inf_state(Mosquito::mState::mInfected);
        }else{
          digesting_mozzies[a]->set_moz_inf_state(Mosquito::mState::mSusceptible);
        }
        //add the mozzie to the current pop and remove from the digeting ones:
        current_mosquitoes.push_back(digesting_mozzies[a]);
        digesting_mozzies.erase(digesting_mozzies.begin()+a);
      }
    }

    //remove the dead mozzies and store the female mozzies who are not dead 
    //2) mozzie dies 
    int nDeaths=0;
    for(int b=0;b<current_mosquitoes.size();b++){
      NumericVector r4=runif(1,0,1);
      if(r4[0]<=daily_pr_death[i]){//
        //set the mozzie state to death 
        current_mosquitoes[b]->death_moz(); 
        ++nDeaths;
 
      }
    }
    
    //remove the dead mozzies 
    for(int k=0; k<current_mosquitoes.size();k++){
      if(current_mosquitoes[k]->get_moz_inf_state()==Mosquito::mState::mDead){
        //remove the dead mozzie:
        current_mosquitoes.erase(current_mosquitoes.begin()+k); 
        //  Rcpp::Rcout << "dead mozzie: " << k << std::endl;
      } 
    }
    
    // Rcpp::Rcout << "female_mozzie_ids: " << female_mozzie_ids << std::endl;
    //0) birth of mozzies:
    int v_new_inf_moz;
    int o_mozzie_births;
      t_rate=mozzie_birth_rate_t[i]; 
      o_mozzie_births=R::rpois(t_rate*dt);
      //add more births as vertical transmitted mozzies:
      if(seed_mozzies==TRUE){
        //generate a random number for mozzies:
        v_new_inf_moz=R::rpois(seed_on[i]*dt); 
      } else{
        v_new_inf_moz=0; 
      }
      mozzie_births=v_new_inf_moz+o_mozzie_births; 
      //total births are equivalant to the number of seeds and the number of eggs
      //(that WILL become eggs)at the the moment. 
      //
      //add the vertical tranmssited mozzies as having no parents,
      //because they are infected anyway and not dependent on the 
      //parents' infection state
      Mosquito* new_mozzie;
      for(int a=0; a<v_new_inf_moz;a++){
        new_mozzie= new Mosquito(Mosquito::mState::mInfected,i,prob_ct_model_change); //create a female mozzie with no parent. 
        //weird since this is vertical transmission!! 
        current_mosquitoes.push_back(new_mozzie);
        new_inf_mozzies++; 
      }
      
        //add new births  as egg states to the population 
        int n=0; 
        //add eggs to a vector separtely and don't add them to the total population yet:
        if(o_mozzie_births>0){
          
          IntegerVector pool = seq_len(current_mosquitoes.size());
          IntegerVector parent_indexx=Rcpp::sample(pool,o_mozzie_births); 
          
          while(n<o_mozzie_births){
            //for(int n=0;n<mozzie_births;n++){
            //for each mozzie that's going to be born, find the female parent:
            int parent_index=parent_indexx[n];
            // int parent_id=female_mozzie_ids[parent_index];
            int parent_id=current_mosquitoes[parent_index]->get_mozzie_id(); 
            //create a new mozzie:
            Mosquito* birth_mozzie;
            Mosquito::mSpecies parent_sp=current_mosquitoes[parent_index]->get_mozzie_species();
            // Rcpp::Rcout << "parent_sp: " << parent_sp << std::endl;
            //get the parent's infection state 
            Mosquito::mState parent_inf= current_mosquitoes[parent_index]->get_moz_inf_state(); 
            birth_mozzie=new Mosquito(parent_id,parent_sp,parent_inf,vertical_transms_prob,i,prob_ct_model_change);//create a mozzie with a parent and their species and their inf state taken to consideration 
            //new mozzie is now in the egg state. set the egg state perios 
            birth_mozzie->set_egg_period(i,egg_to_adult_days); 
            egg_mozzies.push_back(birth_mozzie); 
           // current_mosquitoes.push_back(birth_mozzie);
            //if the birth mozzie is infected, add to the new mozzie counter 
            if(birth_mozzie->get_moz_inf_state()==Mosquito::mState::mInfected){
              new_inf_mozzies++; 
            }
            n++; 
          }
        }
  //  }
    
    //check the egg age, if egg_period>=current time, change them to adult susceptible mozzies, 
    //update the mozzie population 
    //new adults counter:
    int new_adults=0; //start he counter 
    for(int m=0; m<egg_mozzies.size();m++){
      double egg_d =egg_mozzies[m]->get_egg_time();
      if(egg_d>=i){
        egg_mozzies[m]->adult_moz(i, Mosquito::mState::mSusceptible); //make the mozzie a susceptible mozzie 
        //first check if the mozzie should go straight to the overwintering state
        current_mosquitoes.push_back(egg_mozzies[m]); //add the mozzie to the current population 
        //remove the egg from the vector:
        egg_mozzies.erase(egg_mozzies.begin()+m); 
        //add the new adults counter 
        new_adults++; 
      }
    }
    
    //check if the current time is overwintering, if yes send the mozzies to hybernate 
    if(set_the_overwintering_period[i]==1){
      for(int a=0;a<current_mosquitoes.size();a++){
        //if it's winter, then create a random probability send the mozzie to hyberante:
       NumericVector rh=Rcpp::runif(1,0,1);
       if(rh[0]<mozzie_owerwintering_prob[i]){
          current_mosquitoes[a]->mozzie_overwinter(); 
         //set the overwinter age to zero:
         current_mosquitoes[a]->set_overwinter_age(0); 
          //add the mozzie to the hybernating mozzies vector 
          ovewintering_mozzies.push_back(current_mosquitoes[a]); 
          //remove the mozzie from the current population, 
          current_mosquitoes.erase(current_mosquitoes.begin()+a); 
        }
      }
    }
    
    
    //if the time is for the mozzies to become active again, 
    if(active_period[i]==1){
      for(int b=0;b<ovewintering_mozzies.size();b++){
        //generate a random prob:
        NumericVector r_active=Rcpp::runif(1,0,1);
        if(r_active[0]<prob_mozzies_become_active[i]){
          //make the mozzie active and add it to the current pop
          ovewintering_mozzies[b]->mozzie_become_active(); 
          current_mosquitoes.push_back(ovewintering_mozzies[b]);
          ovewintering_mozzies.erase(ovewintering_mozzies.begin()+b);
        }
      }
    }
    
    
    //4)host recovers 
    for(int k=0;k<nhosts;k++){
      if(current_hosts[k]->get_host_inf_state()==Host::hState::hInfected){
        double i_ct=current_hosts[k]->get_current_ct_host(i); 
        double i_tsi=current_hosts[k]->get_infection_start_time_host(); 
        if(i_ct==40){//change to the next state if at the lod
          if(i-i_tsi>4){//at least 4 days days have gone??? 
          current_hosts[k]->recover_host(i_ct);
          }
        }
      }
    }
    
    //death of the hosts:
    //add a counter for dead hosts as well
    int  dead_hosts=0; 
    for(int a=0;a<current_hosts.size();a++){
      //create a random number 
      NumericVector rd=runif(1,0,1);
      if(rd[0]<=prob_death_hosts){
        current_hosts[a]->death_host(); 
        dead_hosts++; 
        //erase the host from the population 
        current_hosts.erase(current_hosts.begin()+a); 
      }
    }
    
    //birth of the hosts 
    int n_birth_hosts= R::rpois(host_birth_rate*dt);
    //create susceptible hosts 
    for(int a=0;a<n_birth_hosts;a++){
      new_host =new Host(); //create a mozzie without parents. all these are females.
      current_hosts.push_back(new_host);
    }
    
    //if the mozzies are overwintering, increase the age:
    if(set_the_overwintering_period[i]==1){
    for(int a=0;a<ovewintering_mozzies.size();a++){
      ovewintering_mozzies[a]->increase_overwinter_age(); 
    }
    }
    
    //if the time period to die if overwintering is greater than or equal to death age; kill the mozzie 
    if(set_the_overwintering_period[i]==1){
      for(int a=0;a<ovewintering_mozzies.size();a++){
          if(ovewintering_mozzies[a]->get_overwinter_age()>=ovewintering_mozzies[a]->get_overwinter_death_age()){
            //kill the mozzie:
            ovewintering_mozzies[a]->death_moz();
            //remove the mozzie from the overwintering mozzies
            ovewintering_mozzies.erase(ovewintering_mozzies.begin()+a); 
          }
      }
    }
    
    
    //increase the age of the living mozzies and those not captured: 
    for(int m=0;m<current_mosquitoes.size();m++){
      current_mosquitoes[m]->increase_mozzie_age(); 
    }
    
    
    //calculate the mozzie deaths and update the mozzie pop. also get the number of s and infected mozzies as well
    IntegerVector nMozzieStates(5);
    for(int a=0;a<current_mosquitoes.size();a++){
      nMozzieStates[0] +=current_mosquitoes[a]->check_moz_sus()*1; //daily prevalence 
      nMozzieStates[1] +=current_mosquitoes[a]->check_moz_inf()*1; //ditto 
    }
  
    nMozzieStates[2]=  nDeaths;//daily deaths 
    nMozzieStates[3] =mozzie_births; //daily birthds 
    nMozzieStates[4]=new_adults;
    timeStateMozzies(i,_)=nMozzieStates;
    
    
    //include the current state to the mmatricses for hosts and mozzies 
    IntegerVector nHostStates(5);
    
    for(int j=0;j<nhosts;j++){
      nHostStates[0] +=current_hosts[j]->check_host_sus()*1;
      nHostStates[1] +=current_hosts[j]->check_host_inf()*1;
      nHostStates[2] +=current_hosts[j]->check_host_rec()*1; 
    }
    
    nHostStates[3]=dead_hosts;
    nHostStates[4]=n_birth_hosts;
    
    //Rcpp::Rcout << "Number of infection hosts: " <<   nHostStates[1] << std::endl;
    timeStatehosts(i,_)=nHostStates;
    
    
    total_mozzies.push_back(current_mosquitoes.size()); 
    new_infs.push_back(new_inf_mozzies); 
    
    //save the ct parameters of the infected mozzes:
    for(int k=0;k<current_mosquitoes.size();k++){
      if(current_mosquitoes[k]->get_moz_inf_state()==1){
        NumericVector moz_pars=current_mosquitoes[k]->get_mozzie_ct_params();
        mozzie_ct_params<< current_mosquitoes[k]->get_mozzie_id() << "," << 
          i << "," <<moz_pars[0] << "," << moz_pars[1] << "," <<
          moz_pars[2] << "," << moz_pars[3] << "," << moz_pars[4] << "\n"; 
      }
    }
    
    //ditto for the hosts:
    for(int k=0;k<current_hosts.size();k++){
      if(current_hosts[k]->get_host_inf_state()==Host::hState::hInfected){
        NumericVector host_pars=current_hosts[k]->get_host_ct_params();
        host_ct_params<<current_hosts[k]->get_host_id() << "," << 
          i << "," << host_pars[0] << "," << host_pars[1] << "," <<
            host_pars[2] << "," << host_pars[3] << "," <<
              host_pars[4] << "\n";
      }
    }
   
    //update the state file in the outputs folders 
    host_states << nHostStates[0]  << "," << nHostStates[1]  << "," << nHostStates[2] 
                << "," << nHostStates[3] << "," << nHostStates[4]  <<  "," << current_hosts.size() << "," << new_inf_hosts << "\n";
    //  mozzie_states<<"Susceptible, Infectious,Deaths, Births, new_adults, Total_pop,\n" ; 
    mozzie_states << nMozzieStates[0] << "," << nMozzieStates[1]  << "," 
                  << nMozzieStates[2]  << "," << nMozzieStates[3]  <<  "," << nMozzieStates[4]  <<  "," <<current_mosquitoes.size()  << ","<< new_inf_mozzies << "\n";
    

    //get the values needed to store the dataset; 
    for(int m=0;m<current_mosquitoes.size();m++){
      double ct_bird_moz=current_mosquitoes[m]->get_bitten_bird_ct_to_mozzie(); 
      double mozzie_ct=40; 
      int ct_type=0; 
      //get the current ct value of the mozzie based on the method:
      if(mozzie_viral_dynamics_method==1){
        ct_type=1;
        mozzie_ct=current_mosquitoes[m]->m1_get_current_viral_load(i);
      }else{
        int ct_ch =current_mosquitoes[m]->get_ct_model_change(); 
        NumericVector mozzie_ct_dat; 
        if(ct_ch==1){
           mozzie_ct_dat=current_mosquitoes[m]->m22_get_current_viral_load(ct_bird_moz);
          
        }else{
           mozzie_ct_dat=current_mosquitoes[m]->m33_get_current_viral_load(ct_bird_moz,prob_ct_model_change,i); 
  
        }

        mozzie_ct=mozzie_ct_dat[0];
        ct_type=mozzie_ct_dat[1];
       
      }
      mozzie_objects_i << current_mosquitoes[m]->get_mozzie_id() << "," <<
        current_mosquitoes[m]->get_moz_inf_state() << "," <<
          i << "," <<   current_mosquitoes[m]->get_inf_start_time() << "," <<
            mozzie_ct << "," <<   ct_type << "," <<
              current_mosquitoes[m]->get_mozzie_bite_counter() << "," <<
                current_mosquitoes[m]->get_mozzie_species() << ","  <<
                  current_mosquitoes[m]->get_mozzie_sex() << "," <<
                    current_mosquitoes[m]->get_mozzie_age() << "," <<
                    current_mosquitoes[m]->get_mozzie_death_age() << "," <<
                sim_number << "\n"; 
    }
    
    mozzie_objects_i.close(); 
    
    //update the hosts states:
    for(int m=0;m<current_hosts.size();m++){
      host_objects_i << current_hosts[m]->get_host_id() << "," <<
        current_hosts[m]->get_host_inf_state() << "," <<
          i << "," <<
            current_hosts[m]->get_infection_start_time_host() << "," <<
              current_hosts[m]->get_current_ct_host(i) << "," <<
              sim_number << "\n";
    }
    host_objects_i.close(); 
  
    i++;
  }
  host_states.close(); 
  mozzie_states.close(); 
  mozzie_ct_params.close(); 
  host_ct_params.close(); 
 
  
}




