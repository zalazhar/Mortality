
get_parameters_CSV <- function(type="age",asis = TRUE) {
  
  parameters <- switch(type,
                       age = read.csv("age_parameters.csv" , as.is = asis),
                       time = read.csv("data/time_parameters.csv" , as.is = asis),
                       correlation = read.csv("data/correlation_matrix.csv" , as.is = asis))
  return (parameters)
  
}


get_parameters <- function(type="age",model = 2014, asis = TRUE) {
  
  # the parameter objects are part of the package and represented as data.frames
  parameters <- get(paste0(type,"_" ,as.character(model)))
  return (parameters)
  
}


noise <- function (num_sims, proj_years, model=2014, seed = 123 ){
 
  corr <- get_parameters(type="correlation", model = model)
  eigen_values <- eigen(corr)
  stoch_variables_names <- c("epsilon_male", 
                             "epsilon_female",
                             "delta_male", 
                             "delta_female")
  
  set.seed(seed)
  norm_simulations <- matrix(rnorm(num_sims*proj_years*length(stoch_variables_names)), 
                             num_sims*proj_years,
                             length(stoch_variables_names))
  simulations <- t(eigen_values$vectors %*% diag(sqrt(eigen_values$values)) %*%t(norm_simulations))
  
  stochastic_variables <- list()
  for (var in 1:length(stoch_variables_names)){
      stochastic_variables[[var]] <- matrix(simulations[,var], num_sims, proj_years)   
  }
  
  return(stochastic_variables)
  
}

time_dynamics <- function (num_sims, proj_years, model = 2014){
  
  types <- list("male.EU" = 1, "female.EU" = 2, "male.NL" = 3, "female.NL" = 4) 
  
  dfTimeParams <- get_parameters("time", model = model)
  
  simulations <- noise(num_sims,proj_years)
  
  dynamic_variables <- list()
  for (var in seq_along(simulations)){
    if (var == types$male.EU  || var == types$female.EU){
      init <- rep(dfTimeParams[var,"K"], max(num_sims,1))
      theta <- dfTimeParams[var,"theta"]
      dynamic_variables[[var]] <- t(apply(cbind(init,theta + simulations[[var]]),1,cumsum))
    }else{
      init <- rep(dfTimeParams[var - 2,"kappa"], max(num_sims,1))
      a <- dfTimeParams[var - 2, "a"]
      dynamic_variables[[var]] <- t(apply(cbind(init, simulations[[var]]),1,
                                          function(x){filter(x,a,"recursive")}))
    }
  }
  
  dynamic_variables <- lapply (dynamic_variables, 
                                function (x) {
                                  colnames (x) <- paste0("y", (model-1):(model + proj_years - 1))
                                  return(x)
                                }
                              )
  return(dynamic_variables)
  
}



qx_lower_90 <- function (num_sims = 1, model = 2014, proj_years = 170){
  
  male_female <- 1:2
  dynamics <- time_dynamics(num_sims,proj_years,model)
  age_parameters <- get_parameters("age", model)
  qx <-list()
  ages <- 0:90
  for (gender in male_female) {
    qx[[gender]] <- lapply(ages, 
                           function(x){1 - exp(-exp(apply_age_factor(age_parameters, x, gender,dynamics)))}
                    )
  }
  
  return(qx)
}

generate_qx <- function (num_sims = 1, model = 2014, proj_years = 170){
  
  qx <- qx_lower_90(num_sims, model, proj_years)
  male_female <- 1:2
  for (gender in male_female){
    total <- 0
    total_weighted <- 0
    for (age in 80:90){
      total <- total - log(-1/log(1-qx[[gender]][[age+1]]) -1)
      total_weighted <- total_weighted - age*log(-1/log(1-qx[[gender]][[age+1]]) -1)
    } 
    for (age in 91:120 ){
      qx[[gender]][[age+1]] <- 
        1 - exp(-1/(1 + exp(-((1/11 - 85*(age - 85)/110)*total + (age - 85)/110*total_weighted))))
    }
    
    #qx[[gender]] <- data.frame(do.call("rbind", qx[[gender]]))
    
    # # convert the list of ages to one big data frame
    # ages <- rep(0:120, each = num_sims)
    # sims <- rep(1:num_sims, 121)
    # qx[[gender]] <- data.frame(do.call("rbind", qx[[gender]]))
    # qx[[gender]]$age <- ages
    # qx[[gender]]$simulation <- sims
    # 
    # #reshape from wide (proection years) to long
    # time_columns <-  names(qx[[gender]])[grep("y",names(qx[[gender]]))]
    # times <- as.numeric(gsub("y","", time_columns))
    # 
    # cat("reshaping ....")
    # qx[[gender]] <- reshape(qx[[gender]], 
    #                         direction = "long", 
    #                         varying = time_columns, 
    #                         v.names = c("value"),
    #                         times = times)
    # 
}
  
  return (qx)
}


get_qx_series <- function(qx, age, gender = "male", start_year = 2014, type ="diagonal"){
  
  qx <- switch(gender,
              male = qx[[1]],
              female = qx[[2]]
        )
  
  projection_years <- start_year - 2 + grep("y", colnames(qx[[1]]))
  end_projection <- max(projection_years)
  end_projection_year <-start_year + end_projection - 1

  qx_series <- list()
  start_age = age
  for (t in start_year:end_projection){
    age <- start_age + t - start_year
    qx_series[[paste0(as.character(t),"_", as.character(age))]] <-
      qx[[min(age,121)]][, paste0("y", as.character(t))]
  }
  
  return(qx_series)
}


apply_age_factor <- function(age_parameters,age, gender,dynamics){
  
  male = 1
  index = age + 1
  if (gender == male){
    age_parameters[index, "male_Ax"] + 
    age_parameters[index, "male_Bx"]*dynamics[[1]]+
    age_parameters[index, "male_alphax"] + 
    age_parameters[index, "male_betax"]*dynamics[[3]]
  }
  else{
    age_parameters[index, "female_Ax"] + 
    age_parameters[index, "female_Bx"]*dynamics[[2]]+
    age_parameters[index, "female_.alphax"] + 
    age_parameters[index, "female_betax"]*dynamics[[4]]
    
  }
  
}

number_of_lifes_series <- function(qx_series){
  start <- 1
  lx <- list()
  lx[[start]] <- 1-qx_series[[start]]
  for (proj in ((start + 1):length(qx_series))){
    lx[[proj]] <- (1-qx_series[[proj]])*lx[[proj - 1]]
  }
  
  return(lx)
}

life_expectancy_series <- function(number_of_lifes_series){
  ex <- list()
  end <- length(number_of_lifes_series)
  ex[[end]] <- number_of_lifes_series[[end]]
  for (proj in (end - 1):1 ){
    
    ex[[proj]] <- number_of_lifes_series[[proj]] + ex[[proj + 1]]
  }
  
  return(ex)
}


# list of projected life expected 
eex_series <- function (qx_series){

    end <- length(qx_series) 
    npx <- list()
    npx[[end]] <-  1- qx_series[[end]] 
    eex <- 0
    
    for (proj in  (end - 1):1){
      npx[[proj]] <- 1-qx_series[[proj]] + npx[[proj + 1]]*(1- qx_series[[proj]])
      eex <- eex + npx[[proj]]
      print (mean(npx))
    }
    
    return(eex)    
  
}

#life expectancy
eex <- function(qx, age, start_year = 2014, gender = "male"){
  
 
  if (gender == "female") gender_index = 2 else gender_index = 1
  
  qx <- qx[[gender_index]]
  number_projections <- length(colnames(qx[[1]]))
  end_projection <- as.numeric(substr(colnames(qx[[1]])[number_projections],2,5))
  
  sumEex <- 0    
  for (t in start_year:end_projection) {
    npx <- 1
    for (projection in start_year:t ){
      parAge <- min(age + projection - start_year +1,120)
      year <- paste0("y", projection)
      npx <- npx*(1-qx[[parAge]][,year])     
    }
   sumEex <- sumEex + npx 
  }
 
  return (sumEex + 0.5)
  
}

ages <- function(){
  return (list(0:120))
}

histLifeExpectancy <- function(eex){
  
  df    <- data.frame(life_expectancy = eex)
  d <- ggplot(df, aes(x = life_expectancy)) +
  geom_histogram(aes(y=..density..),
                colour = "blue", fill = 'white') + 
  geom_density() + geom_vline(aes(xintercept=mean(life_expectancy, na.rm=T)),   # Ignore NA values for mean
                              color="red", linetype="dashed", size=1)   
  d <- d + theme_bw()
  
  return(d)
  
}

drawLifeExpectancyOverTime <- function(startYear,endYear,age){
  
  ls_output <- lapply(startYear:endYear, function(x){mean(eex(preSim,age,x))})
  d <- ggplot(data=data.frame(life_expectancy=unlist(ls_output)), 
              aes(x= c(startYear:endYear), 
              y=life_expectancy, group=1)) +
      geom_line()+
      geom_point()
  
  d<- d+theme_bw()
  
  return(d)
  
}

calculate_quantile <- function(stochast){
  return (data.frame(quantile = round(quantile(stochast, 
                                               probs = c(0.75, 0.80, 0.85, 0.90, 0.95, 0.995)),
                                      2),
                     best_estimate = round(mean(stochast),2)))
}


loadSimData <- function(){
  
  return(preSim)
}


