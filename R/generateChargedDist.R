generateChargedDist <- function(
  proteinSequence = "SEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOW",
  charge = 10, 
  removeFirstAA = FALSE, 
  modification = list( C=0, H=0, N=0, O=0, S=0, P=0, Na=0, K=0, Ca=0)
)
{
  if (removeFirstAA) proteinSequence <- substring(proteinSequence, 2)
  
  proteinChemForm <- pepToForm(proteinSequence)
  
  keys <- unique(c(names(proteinChemForm), names(modification)))
  
  proteinChemForm <- mapply(sum, proteinChemForm[keys], modification[keys])
  proteinChemForm <- setNames(proteinChemForm, keys)
  proteinChemForm <- as.list(proteinChemForm)
  
  dframe <- ldply(charge, function(z, cform){
    
    cform$H <-  cform$H + z
    
    dist <- IsotopicDistribution( cform, charge = z )
    dist$charge <- z
    dist$chemForm <- listToformularChar(cform)  
    return(dist)
    
  }, proteinChemForm)
  
  return(dframe)
}