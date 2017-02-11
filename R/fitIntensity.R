fitIntensity <- function(measuredSpectrum, simulatedSpectrum, mzWidth = 0.2){
  
  simulatedSpectrum_fitted <- ddply(simulatedSpectrum, "charge", function(simSpec, mesSpec){
    
    simSpec100 <- simSpec[ simSpec$percent == 100,]
    filter <- mesSpec$mz > simSpec100$mz - mzWidth/2 &  mesSpec$mz < simSpec100$mz + mzWidth/2
    
    if (any(filter)){
      
      mesSpec <- mesSpec[filter,]
      mesSpec <- mesSpec[ mesSpec$intensity == max(mesSpec$intensity),]
      
      simSpec$intensity <- mesSpec$intensity * simSpec$percent / 100
    }
    else simSpec$intensity <- 0
    
    return(simSpec)
  }, measuredSpectrum)
  
  return(simulatedSpectrum_fitted)
}

