//
// Age structured population model with multiple fleets
//
// Selects R0 and the apical fishing mortality by quarter for each fleet to explain catches and relative depletion. 
// How does it treat MSY...?
//

#include <tuple>
#include <TMB.hpp>
// #include <algorithm>    // std::find
// #include <vector>  

template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
  
}

// see https://stackoverflow.com/questions/35098211/how-do-i-return-two-values-from-a-function

std::tuple<float, float> returns_two() 
{
  return std::make_tuple(1.01, -1.46);
}


// Dan - I'm not sure exactly how this popmodel() within the population model works... 
// I'm guessing this is what you had to do in order to get the thing to return two values, but this is beyond my troubleshooting capabilities for C++
// I'm going to have to assume that there aren't any bugs in the population model itself. 

template<class Type>
std::tuple<Type, Type, Type>  popmodel(Type alpha, 
              Type beta, 
              Type r0, 
              Type ssb,
              int quarters,
              int max_age,
              int recruit_quarters,
              int years,
              vector<Type> weight_at_age,
              vector<Type> maturity_at_age,
              int num_fleets,
              matrix<Type> f_fleet_apical,
              matrix<Type> selectivities,
              matrix<Type> spawning_season_n,
              Type total_ssb0,
              vector<Type> m,
              vector<Type> recruit_data,
              int data_start_quarter,
              vector<Type> depletion,
              matrix<Type> catches
){
  
  matrix<Type> test2(quarters, num_fleets);
  
  Type test = ssb / (alpha+(beta*ssb));
  
  Type tpen = 0.0;
  // number at age
  matrix<Type> n_at_a(quarters,max_age);
  n_at_a = n_at_a.setZero();
  
  // biomass at age (tons)
  matrix<Type> b_at_a(quarters,max_age);
  b_at_a = b_at_a.setZero();
  
  // SSB at age (tons)
  matrix<Type> ssb_at_a(quarters,max_age);
  ssb_at_a = ssb_at_a.setZero();
  
  // Vector of recruits in each quarter
  vector<Type> recruits(quarters);
  recruits = recruits.setZero();
  
  // catch by age (tons)
  matrix<Type> c_at_a(quarters,max_age);
  c_at_a = c_at_a.setZero();
  
  // catch by fleet (tons)
  matrix<Type> c_at_fleet(quarters,num_fleets);
  c_at_fleet = c_at_fleet.setZero();
  
  // sum of squares of catch by fleet (vs. observed)
  matrix<Type> c_at_fleet_ss(quarters,num_fleets);
  c_at_fleet_ss = c_at_fleet_ss.setZero();
  
  // Vector of total catches in each quarter
  vector<Type> c_at_t(quarters);
  c_at_t = c_at_t.setZero();
  
  // Vector of catch sum of squares (vs. observed) in each quarter
  vector<Type> c_ss(quarters);
  c_ss = c_ss.setZero();
  
  // Vector of SSB in each quarter
  vector<Type> ssb_at_t(quarters);
  ssb_at_t = ssb_at_t.setZero();
  
  // Vector of depletion in each quarter
  vector<Type> dep_at_t(quarters);
  dep_at_t = dep_at_t.setZero();
  
  // Vector of depletion in each quarter
  vector<Type> dep_at_q(4);
  dep_at_q = dep_at_q.setZero();
  
  // Vector of avg depletion (SSB) in each year
  vector<Type> dep_avg(years);
  dep_avg = dep_avg.setZero();
  
  // Vector of depletion sum of squares (vs. observed) in each year
  vector<Type> dep_ss(years);
  dep_ss = dep_ss.setZero();
  
  // year counter 
  int year_counter = 0; 
  
  // Vector of recruits in each quarter
  vector<Type> recruit_ss(recruit_quarters);
  recruit_ss = recruit_ss.setZero();
  
  // Recruit counter for recruitment quarters
  int recruit_counter = 0;
  
  // Matrix of f by age in each quarter
  matrix<Type> f_age(quarters, max_age);
  f_age = f_age.setZero();
  
  //////////////////////////////////////////
  
  // Set up time step 1 conditions ------------------------------------------------------
  
  n_at_a(0,0) = r0; // R0 is first age class in first time period
  
  b_at_a(0,0) = n_at_a(0,0) * weight_at_age(0);
  
  ssb_at_a(0,0) = n_at_a(0,0) * maturity_at_age(0) * weight_at_age(0);
  
  recruits(0) = r0;
  
  // ftemp - temporary vector of fishing mortalities to get filled in later? 
  vector<Type> ftemp(num_fleets);
  ftemp = ftemp.setZero();
  
  Type total_f_at_last_t = 0;  // Total fishing mortality in the last time step
  
  
  for (int a = 0; a < max_age; a++){
    
    n_at_a(0,a) = spawning_season_n(a); 
    
    b_at_a(0,a) = n_at_a(0,a) * weight_at_age(a);
    
    ssb_at_a(0, a) =  n_at_a(0,a) * maturity_at_age(a) * weight_at_age(a);
    
  }
  
  ssb_at_t(0) = ssb_at_a.row(0).sum(); // SSB in time 0 is sum of SSB by age in time 0
  
  dep_at_t(0) = ssb_at_t(0);
  dep_at_q(0) = ssb_at_t(0);
  
  // Loop over observed quarters ------------------------------------------------------------
  
  for (int t = 1; t < quarters; t++){ // starting in the second time step
    
    // Use this to define when you want recruitment to happen
    // if ((t%4) >= 0){
      
      n_at_a(t,0) = recruits(t-1);
      
      b_at_a(t,0) = n_at_a(t,0) * weight_at_age(0);
      
    // }
    
    for (int a = 1; a < max_age; a++){
      
      for (int g = 0; g < num_fleets; g++){
        
        ftemp(g) = f_fleet_apical(t-1,g) * selectivities(a-1,g); // calculate f_fleet_age one age at a time. gives a vector of f_fleet_age for each age
        
      }
      
      // Calculate total catches
      c_at_a(t-1, a-1) = weight_at_age(a-1) * n_at_a(t-1, a-1) * (ftemp.sum() / (ftemp.sum() + m(a-1))) * (1 - exp(-(m(a-1) + ftemp.sum())));
      
      // Subtract catches from numbers at age in last time step
      n_at_a(t, a) = posfun(n_at_a(t-1, a-1) * exp(-(m(a-1) + ftemp.sum())), Type(1e-3),tpen); // Grow individuals
      
      b_at_a(t, a) = n_at_a(t, a) * weight_at_age(a);
      
      // Convert numbers to SSB
      ssb_at_a(t, a) =  n_at_a(t, a) * maturity_at_age(a) * weight_at_age(a);
      
    } //close age loop
    
    ssb_at_t(t) = ssb_at_a.row(t).sum();
    
    dep_at_t(t) = ssb_at_t(t);
    
    dep_at_q(t % 4) = ssb_at_t(t);
    
    // Use this to control when spawning happens (currently every quarter)
    // if ((t%4) >= 0){

      // recruits(t) = ssb_at_t(t)/(alpha+(beta*ssb_at_t(t)));
      
      if (recruit_data(t) == 0){
        
        recruits(t) = ssb_at_t(t)/(alpha+(beta*ssb_at_t(t)));
        
      } else {
        
        recruits(t) = recruit_data(t);
        
        recruit_counter = recruit_counter + 1;
      }
      
    // }
    
    // Divy up catches by fleet
    
    c_at_t(t-1) = c_at_a.row(t-1).sum(); //total catch in time t
    
    matrix<Type> f_at_age = f_fleet_apical.row(t-1); // actually f by fleet, used in next line
    
    total_f_at_last_t = f_fleet_apical.row(t-1).sum();
    
    matrix<Type> c_at_fleet_at_age(max_age,num_fleets);
    
    for (int g = 0; g < num_fleets; g++){
      
      Type fgear = (f_fleet_apical(t-1,g) * selectivities.col(g)).sum();
      
      for (int a = 0; a < max_age; a++){
        
        Type total_f_at_a = (f_fleet_apical.row(t-1).array() * selectivities.row(a).array()).sum();
        
        f_age(t-1,a) = total_f_at_a;
        
        Type temp = c_at_a(t-1,a) * ((f_fleet_apical(t-1, g) * selectivities(a,g)) / (total_f_at_a + 1e-3));
        
        c_at_fleet_at_age(a,g) = temp;
        
        // c_at_fleet_at_age(a,g) = posfun(temp,Type(.001),tpen);
        
        
      } // end age loop
      
      c_at_fleet(t-1,g) = c_at_fleet_at_age.col(g).sum();
      
    } // end gear loop
    
   
    // Compute averages in the final year
    
    if(t % 4 == 3){
      
      Type mean = sum(dep_at_q) / 4;
      
      dep_avg(year_counter) = mean;
      
      year_counter = year_counter + 1;
      
      dep_at_q = dep_at_q.setZero();
      
    }
    
  } //close quarters loop
  
  // calculate total catches
  
  vector<Type> quarterly_catch(4);
  
  for (int i = 0 ; i < 3; i++){
    
    quarterly_catch(i) = c_at_fleet.row(quarters - (i + 1)).sum();
   
 }
  
  Type final_annual_catch = quarterly_catch.sum();
  
  // std::cout << final_annual_catch << "\n";
  
  Type final_catch = c_at_fleet.row(quarters - 2).sum();
  
  Type final_ssb = ssb_at_a.row(quarters - 2).sum();
  
  
  // std::cout << final_catch << "\n";
  
  return std::make_tuple(final_catch, final_ssb,final_annual_catch);
}



template<class Type>
Type objective_function<Type>::operator() ()
  {
  
/////////////////
// Load Inputs //
/////////////////
  
  //Load Data
  DATA_INTEGER(max_age);
  DATA_INTEGER(quarters);
  DATA_INTEGER(years);
  DATA_INTEGER(num_fleets);
  DATA_SCALAR(steepness);
  DATA_INTEGER(data_start_quarter);
  DATA_VECTOR(m);
  DATA_VECTOR(s);
  DATA_VECTOR(weight_at_age);
  DATA_VECTOR(maturity_at_age);
  DATA_MATRIX(selectivities);
  DATA_MATRIX(catches);
  DATA_MATRIX(efforts);
  DATA_VECTOR(recruit_data);
  DATA_VECTOR(project_recruit_data);
  DATA_INTEGER(recruit_quarters);
  DATA_VECTOR(depletion);
  DATA_VECTOR(dep_weight);
  DATA_SCALAR(msy);
  DATA_VECTOR(fmult);
  DATA_SCALAR(ssb_msy_to_ssb_0);
  
  // Specify Parameters
  PARAMETER(log_r0);
  PARAMETER_MATRIX(q_est);
  PARAMETER(log_burn_f);
  PARAMETER(log_dep_sigma);
  // PARAMETER(log_catch_sigma);
  // PARAMETER(log_recruit_sigma);

  
  // Type a;
  
  // Type b;
  
  // std::tie(a, b) = returns_two();
  
  // a and b are now 1 and -1, no need to work with std::tuple accessors
  // std::cout << "A" << a << "B is" << b << std::endl;
 
  // Process data and parameters

  Type r0 = exp(log_r0);
  
  matrix<Type> q_at_t(quarters, num_fleets);
  q_at_t = q_at_t.setZero();

  matrix<Type> f_fleet_apical(quarters, num_fleets); // f_fleet_apical is f_fleet_apical
  f_fleet_apical = f_fleet_apical.setZero();
  
  // Deal with fs by fleet
  for (int q = 0; q < quarters; q++){

    for (int g = 0; g < num_fleets; g++){

      if (q < data_start_quarter & g == 0){
        
        f_fleet_apical(q,g) = exp(log_burn_f);
        
      } else{
        
        //f_fleet_apical(q,g) = exp(log_f(q,g)); //convert to raw f_fleet_apical, since I can't figure out how to vectorize
        
        f_fleet_apical(q,g) = exp(q_est(q,g)) * efforts(q,g); //q is calculating f here
        q_at_t(q,g) = exp(q_est(q,g));
          
      }
    } // close fleet loop

  } //close quarter loop


  // population to run to equilibirum
  matrix<Type> zero_n_at_a(100,max_age);
  zero_n_at_a = zero_n_at_a.setZero();
  
  // ssb to run to equilibrium
  matrix<Type> zero_ssb_at_a(100,max_age);
  zero_ssb_at_a = zero_ssb_at_a.setZero();

  // set population in first time and age class to r0
  zero_n_at_a(0,0) = r0;

  /////////////////
  // Equilibrium //
  /////////////////
  
  for (int t=1; t<100; t++){ //equilibrium time loop

    // Use this to define when you want recruitment to happen
    if ((t%4) >= 0){

      zero_n_at_a(t,0) = r0;

    }

    for (int a = 1; a < max_age; a++){ //age loop

      zero_n_at_a(t,a) = zero_n_at_a(t-1, a-1) * exp(-(m(a-1))); // Grow individuals
      
      zero_ssb_at_a(t,a) =  zero_n_at_a(t,a) * maturity_at_age(a) * weight_at_age(a); // Convert numbers to SSB
      
    } //close age loop

  } //close equilibrium time loop

  matrix<Type> spawning_season_ssb = zero_ssb_at_a.row(99);

  matrix<Type> spawning_season_n = zero_n_at_a.row(99);

  Type total_ssb0 = spawning_season_ssb.sum();

  // BH Params
  Type alpha = (total_ssb0 / r0) * (1 - steepness) / ( 4 * steepness);
  
  Type beta = (5 * steepness - 1)/(4 * steepness * r0);
  
  Type pen = 0.0;
  
  //////////////////
  // Data Storage //
  //////////////////
  // Make sure everything actually is zero. Go back and try and use .setZero() here instead.
  // Shouldn't matter though since the problem was in the recruitment
  
  // number at age
  matrix<Type> n_at_a(quarters,max_age);
  n_at_a = n_at_a.setZero();

  // biomass at age (tons)
  matrix<Type> b_at_a(quarters,max_age);
  b_at_a = b_at_a.setZero();
  
  // SSB at age (tons)
  matrix<Type> ssb_at_a(quarters,max_age);
  ssb_at_a = ssb_at_a.setZero();
  
  // Vector of recruits in each quarter
  vector<Type> recruits(quarters);
  recruits = recruits.setZero();

  // catch by age (tons)
  matrix<Type> c_at_a(quarters,max_age);
  c_at_a = c_at_a.setZero();

  // catch by fleet (tons)
  matrix<Type> c_at_fleet(quarters,num_fleets);
  c_at_fleet = c_at_fleet.setZero();

  // sum of squares of catch by fleet (vs. observed)
  matrix<Type> c_at_fleet_ss(quarters,num_fleets);
  c_at_fleet_ss = c_at_fleet_ss.setZero();
  
  // Vector of total catches in each quarter
  vector<Type> c_at_t(quarters);
  c_at_t = c_at_t.setZero();
  
  // Vector of catch sum of squares (vs. observed) in each quarter
  vector<Type> c_ss(quarters);
  c_ss = c_ss.setZero();

  // Vector of SSB in each quarter
  vector<Type> ssb_at_t(quarters);
  ssb_at_t = ssb_at_t.setZero();
  
  // Vector of depletion in each quarter
  vector<Type> dep_at_t(quarters);
  dep_at_t = dep_at_t.setZero();
  
  // Vector of depletion in each quarter
  vector<Type> dep_at_q(4);
  dep_at_q = dep_at_q.setZero();

  // Vector of avg depletion (SSB) in each year
  vector<Type> dep_avg(years);
  dep_avg = dep_avg.setZero();

  // Vector of depletion sum of squares (vs. observed) in each year
  vector<Type> dep_ss(years);
  dep_ss = dep_ss.setZero();
  
  // year counter 
  int year_counter = 0; 
  
  // Vector of recruits in each quarter
  vector<Type> recruit_ss(recruit_quarters);
  recruit_ss = recruit_ss.setZero();
  
  // Recruit counter for recruitment quarters
  int recruit_counter = 0;
  
  // Matrix of f by age in each quarter
  matrix<Type> f_age(quarters, max_age);
  f_age = f_age.setZero();
  
//////////////////////////////////////////

// Set up time step 1 conditions ------------------------------------------------------
  n_at_a(0,0) = r0; // R0 is first age class in first time period
  
  b_at_a(0,0) = n_at_a(0,0) * weight_at_age(0);
  
  ssb_at_a(0,0) = n_at_a(0,0) * maturity_at_age(0) * weight_at_age(0);

  recruits(0) = r0;

  // ftemp 
  vector<Type> ftemp(num_fleets);
  ftemp = ftemp.setZero();
  
  Type total_f_at_last_t = 0;  // Not sure what this does

  for (int a = 0; a < max_age; a++){

    n_at_a(0,a) = spawning_season_n(a); // numbers at age are the numbers in that spawning season?
    
    b_at_a(0,a) = n_at_a(0,a) * weight_at_age(a);

    ssb_at_a(0, a) =  n_at_a(0,a) * maturity_at_age(a) * weight_at_age(a);
    
  }

  ssb_at_t(0) = ssb_at_a.row(0).sum(); // SSB in time 0 is sum of SSB by age in time 0

  dep_at_t(0) = ssb_at_t(0) / total_ssb0;
  dep_at_q(0) = ssb_at_t(0) / total_ssb0;

// Loop over observed quarters ------------------------------------------------------------

  for (int t = 1; t < quarters; t++){ // starting in the second time step

    // Use this to define when you want recruitment to happen
    if ((t%4) >= 0){

    n_at_a(t,0) = recruits(t-1);
      
    b_at_a(t,0) = n_at_a(t,0) * weight_at_age(0);

    }

    for (int a = 1; a < max_age; a++){

      for (int g = 0; g < num_fleets; g++){

        ftemp(g) = f_fleet_apical(t-1,g) * selectivities(a-1,g); // calculate f_fleet_age one age at a time. gives a vector of f_fleet_age for each age
        
      }
      
      // Calculate total catches
      c_at_a(t-1, a-1) = weight_at_age(a-1) * n_at_a(t-1, a-1) * (ftemp.sum() / (ftemp.sum() + m(a-1))) * (1 - exp(-(m(a-1) + ftemp.sum())));
      
      // Subtract catches from numbers at age in last time step
      n_at_a(t, a) = posfun(n_at_a(t-1, a-1) * exp(-(m(a-1) + ftemp.sum())), Type(1e-3),pen); // Grow individuals
      
      b_at_a(t, a) = n_at_a(t, a) * weight_at_age(a);

      // Convert numbers to SSB
      ssb_at_a(t, a) =  n_at_a(t, a) * maturity_at_age(a) * weight_at_age(a);

    } //close age loop

    ssb_at_t(t) = ssb_at_a.row(t).sum();
    
    dep_at_t(t) = ssb_at_t(t) / total_ssb0;
    dep_at_q(t % 4) = ssb_at_t(t) / total_ssb0;

    // Use this to control when spawning happens (currently every quarter)
    // if ((t%4) >= 0){
      
      recruits(t) = ssb_at_t(t)/(alpha+(beta*ssb_at_t(t)));
      
        if (recruit_data(t) == 0){
          
          recruits(t) = ssb_at_t(t)/(alpha+(beta*ssb_at_t(t)));
          
        } else {
      
          recruits(t) = recruit_data(t);
      
          // Type nll_recruit = -(dnorm(log(recruit_data(t)), log(ssb_at_t(t)/(alpha+(beta*ssb_at_t(t)))), exp(log_recruit_sigma), true));
          // recruit_ss(recruit_counter) = nll_recruit;
          // 
          recruit_counter = recruit_counter + 1;
        }
      
    // }

    // Divy up catches by fleet
    
    c_at_t(t-1) = c_at_a.row(t-1).sum(); //total catch in time t

    matrix<Type> f_at_age = f_fleet_apical.row(t-1); // actually f by fleet, used in next line

    total_f_at_last_t = f_fleet_apical.row(t-1).sum();
    
    matrix<Type> c_at_fleet_at_age(max_age,num_fleets);

    for (int g = 0; g < num_fleets; g++){

      Type fgear = (f_fleet_apical(t-1,g) * selectivities.col(g)).sum();

      for (int a = 0; a < max_age; a++){
        
        Type total_f_at_a = (f_fleet_apical.row(t-1).array() * selectivities.row(a).array()).sum();
        
        f_age(t-1,a) = total_f_at_a;
        
        Type temp = c_at_a(t-1,a) * ((f_fleet_apical(t-1, g) * selectivities(a,g)) / (total_f_at_a + 1e-3));
        
        c_at_fleet_at_age(a,g) = temp; 
        
        // c_at_fleet_at_age(a,g) = posfun(temp, Type(0.001), penalty);
        
      } // end age loop

      c_at_fleet(t-1,g) = c_at_fleet_at_age.col(g).sum();
      
      if (t > data_start_quarter) {
        
      // std::cout <<  c_at_fleet(t-1,g) << "\n";
        
      // Type nll_catch = -(dnorm(log(catches(t-1-data_start_quarter, g) + 1e-3), log(c_at_fleet(t-1, g)), exp(log_catch_sigma), true));
     
     Type nll_catch = -(dnorm(log(catches(t-1-data_start_quarter, g) + 1e-3), log(c_at_fleet(t-1, g) + 1e-3), Type(0.01), true));
     
      c_at_fleet_ss(t-1,g) = nll_catch;
      
      }
    } // end gear loop

  // Summarize time step by computing catch and depletion ss
  
    if (t > data_start_quarter) {
      
      c_ss(t-1) = c_at_fleet_ss.row(t-1).sum(); //catch sum of squares
      
    }
    
    //dep_ss(t) = pow(log(dep_at_t(t) + .01) - log(depletion(t) + 01),2); // calculate depletion sum of squares

    // Compute averages in the final year
    
    if(t % 4 == 3){
      
      Type mean = sum(dep_at_q) / 4;
      
      dep_avg(year_counter) = posfun(mean, Type(1e-3),pen);
      
      dep_ss(year_counter) = pow(log(dep_avg(year_counter) + .01) - log(depletion(year_counter) + .01), 2);
      
      year_counter = year_counter + 1;
      
      dep_at_q = dep_at_q.setZero();
      
    }
    
  } //close quarters loop

// Fill in the final year's catch. There has to be a better way ----------------------------------------------------

int t = quarters - 1; // Account for zero indec

for (int a = 0; a < max_age; a++)
{
  for (int g = 0; g < num_fleets; g++){

    ftemp(g) = f_fleet_apical(t,g) * selectivities(a,g); // calculate f_fleet_apical at age for gear in quarter

    //std::cout << ftemp(g) << "\n";
  }
  
//Calculate catch at age in time t

  c_at_a(t,a) = weight_at_age(a) * (ftemp.sum() / (ftemp.sum() + m(a))) *  n_at_a(t, a ) * (1 - exp(-(m(a) + ftemp.sum())));
}

c_at_t(t) = c_at_a.row(t).sum(); // total catch in time t

matrix<Type> c_at_fleet_at_age(max_age,num_fleets);

// Divide catch by fleet in each time period
for (int g = 0; g < num_fleets; g++){

  Type fgear = (f_fleet_apical(t,g) * selectivities.col(g)).sum();

  for (int a = 0; a < max_age; a++){

    Type total_f_at_a = (f_fleet_apical.row(t).array() * selectivities.row(a).array()).sum();
    
    // if (total_f_at_a <= 1e-3){
    //  std::cout << total_f_at_a << "\n";
    //   
    // }
    
      f_age(t,a) = total_f_at_a;
      
      Type temp = c_at_a(t-1,a) * ((f_fleet_apical(t-1, g) * selectivities(a,g)) / (total_f_at_a + 1e-3));
      
      c_at_fleet_at_age(a,g) = temp;
      
      // c_at_fleet_at_age(a,g) = posfun(temp, Type(1e-3), penalty);
      
      // if ( c_at_fleet_at_age(a,g) <= 10){
      //  
      //  std::cout << c_at_fleet_at_age(a,g) << "\n";
      // 
      // }
      
  }

  c_at_fleet(t,g) = c_at_fleet_at_age.col(g).sum();

  // if (  c_at_fleet(t,g) <=0 ){
    
    // std::cout <<  c_at_fleet(t,g) << "\n";
    
  // }
  
  Type nll_catch = -(dnorm(log(catches(t - data_start_quarter, g) + 1e-3), log(c_at_fleet(t, g) + 1e-3), Type(0.01), true));
  
  // Type nll_catch = -(dnorm(log(catches(t - data_start_quarter, g) + 1e-3), log(c_at_fleet(t, g)), exp(log_catch_sigma), true));
  
  c_at_fleet_ss(t,g) = nll_catch;

}

c_ss(t) = c_at_fleet_ss.row(t).sum(); //catch sum of squares

Type b0 = b_at_a.row(0).sum();

Type n0 = n_at_a.row(0).sum();

Type rtest = recruits(0);

// calculate MSY

vector<Type> mean_recent_f(num_fleets); // eventually make this mean over the last few years, but for just do last year

matrix<Type> recent_f(4 * 5,num_fleets);

recent_f = f_fleet_apical.block(quarters - 4*5 + -1,0,4 * 5,num_fleets);

// std::cout << recent_f << "\n";

for (int g = 0; g < num_fleets; g++){

  // mean_recent_f(g) = f_fleet_apical(quarters - 1, g);
  
  mean_recent_f(g) = recent_f.col(g).mean();
    
}

matrix<Type> f_fleet_recent(quarters, num_fleets); // f_fleet_apical is f_fleet_apical

f_fleet_recent = f_fleet_recent.setZero();

// fill in base case f
for (int i = 0; i < quarters; i++){
  
  for (int g = 0; g < num_fleets; g ++){
    
    f_fleet_recent(i,g) = mean_recent_f(g);
    
  }
  
}

vector<Type> fmult_catches(fmult.size());

vector<Type> fmult_ssb(fmult.size());

// recruit_data =  recruit_data.setZero();

// std::cout << f_fleet_recent << "\n";

for (int i = 0; i < fmult.size(); i++){
 
matrix<Type> f_temp(quarters, num_fleets);

// recruit_data =  recruit_data.setZero();

f_temp = f_fleet_recent * fmult(i);


Type final_catch;
Type final_ssb;
Type final_annual_catch;


std::tie(final_catch,final_ssb,final_annual_catch) = popmodel(alpha,
              beta,
              r0, 
              total_ssb0, 
              quarters, 
              max_age,
              recruit_quarters,
              years,
              weight_at_age, 
              maturity_at_age,
              num_fleets,
              f_temp,
              selectivities,
              spawning_season_n,
              total_ssb0,
              m,
              project_recruit_data,
              data_start_quarter,
              depletion,
              catches);

fmult_catches(i) = final_annual_catch;

fmult_ssb(i) = final_ssb;

 }

// std::cout << "MSY is" << max(fmult_catches) << "\n";
 
 int msy_location = 0;

for (int i = 0; i < fmult_catches.size(); i++){
  if (fmult_catches(i) == max(fmult_catches)){

    msy_location = i;

  }

}

Type msy_hat = max(fmult_catches);

Type msy_penalty = pow(log(msy) - log(msy_hat),2);

vector<Type> fmsy(num_fleets);

fmsy = mean_recent_f * fmult(msy_location);

Type ssb_msy = fmult_ssb(msy_location);

Type ssb_msy_to_ssb_0_hat = ssb_msy / total_ssb0;

Type ssbmsy_penalty = pow(ssb_msy_to_ssb_0_hat - ssb_msy_to_ssb_0,2);


// fmsy = mean_recent_f * fmult;

// std::cout<< "msy location is" << msy_location << "\n";

  // REPORT THINGs!

  
  REPORT(f_fleet_recent);
  
  REPORT(msy_location);
  
  REPORT(fmult_catches);
  
  REPORT(fmult_ssb);
  
  REPORT(ssb_msy);
  
  REPORT(q_at_t);
  
  REPORT(dep_at_t);
  
  REPORT(dep_avg);

  REPORT(dep_ss);

  REPORT(total_ssb0);

  REPORT(alpha);

  REPORT(beta);

  REPORT(r0);
  
  REPORT(rtest);

  REPORT(n_at_a);
  
  REPORT(b_at_a);
  
  REPORT(b0);
  
  REPORT(n0);
  
  REPORT(ssb_at_t);
  
  REPORT(ssb_at_a);
  
  REPORT(recruits);

  REPORT(c_at_fleet); // catch by fleet over time
  
  REPORT(c_at_a); // catch by age over time

  REPORT(c_ss); //sum of squares between fleet catches observed and predicted summed across all fleets

  REPORT(f_fleet_apical); // apical fs by time and fleet
  
  REPORT(f_age); //fs by age over time
  
  REPORT(fmult);
  
  REPORT(fmsy);
  
  
  REPORT(msy_hat);

  // Type ss = sum(dep_ss)*3 + sum(c_ss); // Return total sum of squares, no weighting
  
  Type nll_dep = -sum(dnorm(log(depletion), log(dep_avg), exp(log_dep_sigma), true));

  //Type ll_catch = -sum(dnorm(log(depletion), log(dep_at_t), Type(0.05), true));

  Type total_nll_c = sum(c_ss);
  
  Type fpen = pow(Type(1.57) - sum(mean_recent_f) / sum(fmsy),2);

  Type bpen = pow(Type(1.0) - ssb_at_t(quarters - 1) / ssb_msy,2);
  
  Type total_nll = (dep_weight(0) * nll_dep) + (dep_weight(1) * total_nll_c) + pen;// + msy_penalty + fpen + bpen;
  
  return(total_nll);

}
