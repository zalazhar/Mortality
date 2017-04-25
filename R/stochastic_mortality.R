library(ggplot2)

getParameters_CSV <- function(type="age",asis = TRUE) {
  
  parameters <- switch(type,
                       age = read.csv("age_parameters.csv" , as.is = asis),
                       time = read.csv("data/time_parameters.csv" , as.is = asis),
                       correlation = read.csv("data/correlation_matrix.csv" , as.is = asis))
  return (parameters)
  
}

getParameters <- function(type="age",asis = TRUE) {
  
  parameters <- switch(type,
                       age = age,
                       time = time,
                       correlation = correlation)
  return (parameters)
  
}

kappa <- function (numberSims = 0, endYear = 2184) {
  
  startYear = 2013
  if (endYear <= 2013) {endYear = 2013} 
  # number of projection years
  progYears <- endYear - startYear + 1
  
  noise <- noise(numberSims,progYears)
  
  # the levels in which we need the parameters for creating the mortality intensity
  male = 1 ; female = 2 ;EU = 1;NL = 2
  dfTimeParams <- getParameters("time")
  
  listOutput <- list()
  for (gender in male:female) { 
    for (region in NL:EU) { 
      type <- paste0(gender,region)
      df <- data.frame()
      for (year in 1:progYears) {
          indexYear = paste0("y",year + startYear - 1)
          indexYearPrev = paste0("y",year + startYear - 2)
          if (year == 1){
            startTimeParameters <- rep(dfTimeParams[gender,region + 1], max(numberSims,1))
            df <- data.frame("y2013" = startTimeParameters)} else{
              if (region == EU) {
                df[,indexYear] <-df[,indexYearPrev] + dfTimeParams[gender,"theta"] + noise$epsilon[[gender]][,year]
              }else{
                df[,indexYear] <- dfTimeParams[gender,"a"]*df[,indexYearPrev] + noise$delta[[gender]][,year]      
              }   
           }  
      } 
   listOutput[[type]] <- df
 }

}
return (listOutput)
}

noise  <- function(numSims,progYears){
 
  Z1 <- matrix(0,1,progYears);Z2 <- matrix(0,1,progYears)
  if (numSims != 0) {Z1 <- matrix(rnorm(numSims*progYears),numSims,progYears)
                     Z2 <- matrix(rnorm(numSims*progYears),numSims,progYears)}
  
  male = 1;female = 2
  corr <- getParameters("correlation")
  
  epsilon <- list();cov <- list();delta <- list() 
  
  for (gender in male:female){
    epsilon[[gender]] <- sqrt(corr[gender,"variance_e1"])*Z1
    cov[[gender]]  <- corr[gender,"covariance"]/(sqrt(corr[gender,"variance_e1"])*sqrt(corr[gender,"variance_e2"]))
    delta[[gender]] <-  sqrt(corr[gender,"variance_e2"])*(cov[[gender]]*Z1 + sqrt(1 - cov[[gender]]^2)*Z2 )
  }
  
  return (list(epsilon = epsilon, delta = delta))
  
}


listQx <- function( age = 0, numberSims = 0, progYear = 2184){
  
  # making list of dataframes, whereby each item of the list represents a specific age, and
  # each age (list item) contains the kappa and k simulations for male and female as well.
  # '11' = NL & Male, '12' = NL & Female, '21' = EU & Male, '22' = EU & Female
  if (age > 90) {age = 80}
  maleEU = '11'; maleNL = '12' ; FemaleEU = '21' ; FemaleNL = '22'
  maxAge = 120
  ageParameters <- getParameters("age")
  index = maxAge - age + 1; index90 = 90 - age + 1 ;index80 = 80 - age + 1 
  
  kappaPar <- kappa(numberSims, progYear)
  
  qxMale <- list(); qxFemale <- list(); uxMale <- list();uxFemale <- list()
  for (x in 1:index90) { 
     parAge <- age  + x - 1
     lnUxMaleEU <- ageParameters[parAge + 1,"male_Ax"] + ageParameters[parAge+ 1,"male_Bx"]*kappaPar[[maleEU]]
     lnUxMaleNL <- ageParameters[parAge +1,"male_alphax"] + ageParameters[parAge+1,"male_betax"]*kappaPar[[maleNL]]
     lnUxFemaleEU <- ageParameters[parAge +1,"female_Ax"] + ageParameters[parAge+1,"female_Bx"]*kappaPar[[FemaleEU]]
     lnUxFemaleNL <- ageParameters[parAge+1,"female_.alphax"] + ageParameters[parAge+1,"female_betax"]*kappaPar[[FemaleNL]]
     uxMale[[x]] <- exp(lnUxMaleNL + lnUxMaleEU)
     uxFemale[[x]] <- exp(lnUxFemaleNL + lnUxFemaleEU)
     qxMale[[x]] <- 1 - exp(-uxMale[[x]])
     qxFemale[[x]] <- 1 - exp(-uxFemale[[x]])
     qxMale[[x]]$age <- age + x -1
     qxFemale[[x]]$age <- age + x -1
  }
  
  #for ages > 90 we use an extrapolation method based on the mortality intensities of age-range 80-90
  sumhkMale <- 0 ; sumhkFemale <- 0 ;sumhkYMale <- 0 ; sumhkYFemale <- 0
  for (k in 1:11){
    indexk <- index80 + k - 1
    yk <- age + indexk - 1
    sumhkMale <- sumhkMale - log(1/uxMale[[indexk]] - 1 )
    sumhkFemale <- sumhkFemale - log(1/uxFemale[[indexk]] - 1) 
    sumhkYMale <- sumhkYMale - yk*(log(1/uxMale[[indexk]] - 1 ))
    sumhkYFemale <- sumhkYFemale - yk*(log(1/uxFemale[[indexk]] - 1 ))
  }
  
 
  startIndex <- index90 + 1
  for (x in startIndex:index) {
     parAge <- age + x - 1   
     gxMale <- (1/11 - 85*(parAge - 85)/110)*sumhkMale + (parAge - 85)/110*sumhkYMale
     gxFemale <- (1/11 - 85*(parAge - 85)/110)*sumhkFemale + (parAge - 85)/110*sumhkYFemale
     uxMale[[x]] <- 1/(1 + exp(-gxMale)) ; uxFemale[[x]] <- 1/(1 + exp(-gxFemale))
     qxMale[[x]] <- 1 - exp(-uxMale[[x]]) ; qxFemale[[x]] <- 1 - exp(-uxFemale[[x]])
     qxMale[[x]]$age <- parAge ; qxFemale[[x]]$age <- parAge
  }
  

  return (list(qxMale = qxMale, qxFemale = qxFemale))
  
}

writeQx <- function(listQx){
  write.csv(listQx,"klm.csv")
  return (invisible)
  
}



#actuarial function ex, the life expectancy expressed in years
#input is the stochastic Qx dataset, which can be generated with the function listQx
#this stochastic dataset contains 120 age-items, each item having n simulations for t 
#projection years

eex <- function(listQx, age, startYear = 2014, gender = "Male"){
  
  qx <- switch(gender , Male = listQx$qxMale,Female = listQx$qxFemale)
  yearProg <- 2184
  
  sumEex <- 0    
  for (t in startYear:yearProg) {
    npx <- 1
    for (prog in startYear:t ){
      parAge <- min(age + prog - startYear+1,121)
      year <- paste0("y", prog)
      npx <- npx*(1-qx[[parAge]][,year])     
    }
   sumEex <- sumEex + npx 
  }
 
  return (sumEex + 0.5)
  
}

ages <- function(){
  return (list(0:120))
}

drawLifeExpectancy <- function(eex){
  
  
  df    <- data.frame(life_expectancy = eex)
  d <- ggplot(df, aes(x = life_expectancy)) +
  geom_histogram(aes(y=..density..),
                colour = "blue", fill = 'white') + 
  geom_density() + geom_vline(aes(xintercept=mean(life_expectancy, na.rm=T)),   # Ignore NA values for mean
                              color="red", linetype="dashed", size=1)   
  d <- d + theme_bw()
  
  return(d)
  
}

calculate_quantile <- function(stochast){
  return (data.frame(quantile = quantile(stochast, probs = c(0.75, 0.80, 0.85, 0.90, 0.95, 0.995))))
}
getDataTable <- function(dt){
  return (datatable(df))
}

loadSimData <- function(){
  
  return(preSim)
}

