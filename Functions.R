

pepToForm <- function(sequence, output = "vector"){
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
    if (output == "vector") result <- paste("C",resultsVector[1],"H",resultsVector[2],"N",resultsVector[3],"O",resultsVector[4],"S",resultsVector[5])
    if (output == "list") result <- list(C = resultsVector[1], H = resultsVector[2], N = resultsVector[3], O = resultsVector[4], S = resultsVector[5])
    return(result)
  }
  else return(NA)
}

addVariableModification <- function(){
  
  
}


# This function is based on the function IsotopicDistribution() from the package OrgMassSpecR
# https://github.com/OrgMassSpec/OrgMassSpecR
# I altered it to fit in my code / programming style

isotopicDist <- function (formula = list(), charge = 1) 
{
  if (charge == 0) 
    stop("a charge of zero is not allowed")
  inputFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0, 
                       Br = 0, Cl = 0, F = 0, Si = 0)
  inputFormula[names(formula)] <- formula
  simulation <- function(inputFormula) {
    massCarbon <- sum(sample(c(12, 13.0033548378), size = inputFormula$C, 
                             replace = TRUE, prob = c(0.9893, 0.0107)))
    massHydrogen <- sum(sample(c(1.0078250321, 2.014101778), 
                               size = inputFormula$H, replace = TRUE, prob = c(0.999885, 
                                                                               0.000115)))
    massNitrogen <- sum(sample(c(14.0030740052, 15.0001088984), 
                               size = inputFormula$N, replace = TRUE, prob = c(0.99632, 
                                                                               0.00368)))
    massOxygen <- sum(sample(c(15.9949146221, 16.9991315, 
                               17.9991604), size = inputFormula$O, replace = TRUE, 
                             prob = c(0.99757, 0.00038, 0.00205)))
    massSulfer <- sum(sample(c(31.97207069, 32.9714585, 33.96786683, 
                               35.96708088), size = inputFormula$S, replace = TRUE, 
                             prob = c(0.9493, 0.0076, 0.0429, 2e-04)))
    massPhosphorus <- inputFormula$P * 30.97376151
    massBromine <- sum(sample(c(78.9183376, 80.916291), size = inputFormula$Br, 
                              replace = TRUE, prob = c(0.5069, 0.4931)))
    massChlorine <- sum(sample(c(34.96885271, 36.9659026), 
                               size = inputFormula$Cl, replace = TRUE, prob = c(0.7578, 
                                                                                0.2422)))
    massFluorine <- inputFormula$F * 18.9984032
    massSilicon <- sum(sample(c(27.9769265327, 28.97649472, 
                                29.97377022), size = inputFormula$Si, replace = TRUE, 
                              prob = c(0.922297, 0.046832, 0.030872)))
    massMolecule <- sum(massCarbon, massHydrogen, massNitrogen, 
                        massOxygen, massSulfer, massPhosphorus, massBromine, 
                        massChlorine, massFluorine, massSilicon)
    mz <- massMolecule/abs(charge)
    return(mz)
  }
  sim <- replicate(10000, expr = simulation(inputFormula))
  b <- seq(from = min(sim) - (1/(2 * abs(charge))), to = max(sim) + 
             1, by = 1/abs(charge))
  bins <- cut(sim, breaks = b)
  mz <- round(tapply(sim, bins, mean), digits = 2)
  intensity <- as.vector(table(bins))
  spectrum <- data.frame(mz, intensity)
  spectrum <- spectrum[spectrum$intensity != 0, ]
  spectrum$percent <- with(spectrum, round(intensity/max(intensity) * 
                                             100, digits = 2))
  row.names(spectrum) <- 1:(nrow(spectrum))
  return(spectrum)
}
