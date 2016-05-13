drawLifeExpectancy <-
function(age){
  qxStochasticSet <- readRDS('data/QxStochasticSet.rds')
  df    <- data.frame(life_expectancy = eex(qxStochasticSet,age))
  ggplot(df, aes(x = life_expectancy)) +
  geom_histogram(aes(y=..density..),
                 binwidth= 0.5, colour = "blue", fill = 'white') + 
  geom_density() + geom_vline(aes(xintercept=mean(life_expectancy, na.rm=T)),   # Ignore NA values for mean
                              color="red", linetype="dashed", size=1)   
  
}
