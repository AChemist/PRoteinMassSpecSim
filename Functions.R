

pepToForm <- function(sequence){
  if (!is.na(sequence)) {
    peptideVector <- strsplit(sequence, split = "")[[1]]
    FindElement <- function(residue){
      element <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      if (residue == "A") element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)
      if (residue == "R") element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)
      if (residue == "N") element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)
      if (residue == "D") element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)
      if (residue == "C") element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)
      if (residue == "c") element <- c(C = 5, H = 6, N = 2, O = 2, S = 1)
      if (residue == "E") element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)
      if (residue == "Q") element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)
      if (residue == "G") element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)
      if (residue == "H") element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)
      if (residue == "I") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if (residue == "L") element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)
      if (residue == "K") element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)
      if (residue == "M") element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)
      if (residue == "m") element <- c(C = 5, H = 9, N = 1, O = 2, S = 1)
      if (residue == "F") element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)
      if (residue == "P") element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)
      if (residue == "S") element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)
      if (residue == "T") element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)
      if (residue == "W") element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)
      if (residue == "Y") element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)
      if (residue == "V") element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)
      return(element)
    }
    resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
    for (i in 1:length(peptideVector)) {
      resultsVector <- FindElement(peptideVector[i]) + resultsVector
    }
    resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, O = 1, S = 0)
    
    result <- paste("C",resultsVector[1],"H",resultsVector[2],"N",resultsVector[3],"O",resultsVector[4],"S",resultsVector[5], "P",0)
    
    return(result)
  }
  else return(NA)
}


formulaCharToList <- function(char){
  
  vec <- unlist(strsplit(as.character(char), split = " "))
  lis <- list(
    C = as.integer(vec[2]),
    H = as.integer(vec[4]),
    N = as.integer(vec[6]),
    O = as.integer(vec[8]),
    S = as.integer(vec[10]),
    P = as.integer(vec[12])
  )
  return(lis)
}

listToformularChar <- function(lis){
  formChar <- paste("C",lis[1],"H",lis[2],"N",lis[3],"O",lis[4],"S",lis[5],"P",lis[6])
  return(formChar)
}


formulaCharToNamedVec <- function(char){
  
  vec <- unlist(strsplit(as.character(char), split = " "))
  vec <- c(
    C = as.integer(vec[2]),
    H = as.integer(vec[4]),
    N = as.integer(vec[6]),
    O = as.integer(vec[8]),
    S = as.integer(vec[10]),
    P = as.integer(vec[12])
  )
  return(vec)
}

namedVecToformularChar <- function(nvec){
  formChar <- paste("C",nvec[1],"H",nvec[2],"N",nvec[3],"O",nvec[4],"S",nvec[5],"P",nvec[6])
  return(formChar)
}


fetchProteinSequence <- function(
  uniprotSpeciesName = "Homo sapiens",
  proteinAccession = "P02144"
)
{
  
  tmp <- availableUniprotSpecies()
  
  if (any(tmp$`Species name` == uniprotSpeciesName)) taxId <- as.integer(tmp$`taxon ID`[tmp$`Species name` == uniprotSpeciesName])
  #else generate warning and quit
  
  tmp <- UniProt.ws(taxId)
  
  proteinSequence <- select(tmp, proteinAccession, "SEQUENCE", "UNIPROTKB")[,2]

  return(proteinSequence)
}


generateChargedDist <- function(
  proteinSequence = "SEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOWSEBASTIANMALCHOW",
  charge = 10, 
  removeFirstAA = FALSE, 
  modification = "C 0 H 0 N 0 O 0 S 0 P 0"
)
{
  if (removeFirstAA) proteinSequence <- substring(proteinSequence, 2)
  
  proteinChemForm <- pepToForm(proteinSequence)
  
  proteinChemForm <- formulaCharToNamedVec(proteinChemForm) + formulaCharToNamedVec(modification)
  proteinChemForm <- namedVecToformularChar(proteinChemForm)
  
  dframe <- ldply(charge, function(x){
    
    tmp <- formulaCharToNamedVec(proteinChemForm)
    tmp <- unname(tmp + c(0,x,0,0,0,0))
    tmp <- list( C = tmp[1], H = tmp[2], N = tmp[3], O = tmp[4], S = tmp[5], P = tmp[6])
    dist <- IsotopicDistribution( tmp, charge = x )
    dist$charge <- x
    dist$chemForm <- listToformularChar(tmp)  
    return(dist)
    
  })
  
  return(dframe)
}

fitIntensity <- function(measuredSpectrum, simulatedSpectrum, decimalPlaces = 2){
  
  fit <- measuredSpectrum[ round(measuredSpectrum$mz, digits = decimalPlaces) %in% round(simulatedSpectrum$mz[ simulatedSpectrum$percent == 100 ], digits = decimalPlaces),]
  fit <- fit[ fit$intensity == max(fit$intensity),]
  
  simulatedSpectrum$intensity <- fit$intensity * simulatedSpectrum$percent / 100  
  
  return(simulatedSpectrum)
  
}
  
