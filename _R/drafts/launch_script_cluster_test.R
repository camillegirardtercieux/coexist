launch_model <- function(mortality, fecundity, seed, seed_r){
  print(mortality)
  print(fecundity)
  print(seed)
  print(seed_r)
  
  if (mortality =="fixed"){
    mortality_fixed<-TRUE
    mortality_proportion<-FALSE
    mortality_stocha<-FALSE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="prop"){
    mortality_fixed<-FALSE
    mortality_proportion<-TRUE
    mortality_stocha<-FALSE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="stocha"){
    mortality_fixed<-FALSE
    mortality_proportion<-FALSE
    mortality_stocha<-TRUE
    mortality_stocha_basal<-FALSE
  }
  
  if (mortality =="stocha_basal"){
    mortality_fixed<-FALSE
    mortality_proportion<-FALSE
    mortality_stocha<-FALSE
    mortality_stocha_basal<-TRUE
  }
  
  if(fecundity=="fixed"){
    nb_seeds_dep_abund<-FALSE
    nb_seeds_indep_abund<-TRUE
  }
  
  if(fecundity=="abund"){
    nb_seeds_dep_abund<-TRUE
    nb_seeds_indep_abund<-FALSE
  }
  
  print(mortality_stocha_basal)
  print(mortality_stocha)
  print(mortality_proportion)
  print(mortality_fixed)
  print(nb_seeds_dep_abund)
  print(nb_seeds_indep_abund)
  
  source(file = here::here("_R", "Basic_parameters_cluster.R"), local=TRUE)
  
 
  print(mod)
  print(nb_seeds)
  print(mort)
  print(start)
  print(model)
  print(model_perf)

}