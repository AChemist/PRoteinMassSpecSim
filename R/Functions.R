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

lis <- proteinChemForm

listToformularChar <- function(lis){
  
  vec <- unlist(lis)
  formChar <- ""
  
  for (i in length(vec)){
    if (vec[i] != 0)  formChar <- paste0(formChar, names(vec)[i], vec[i], " ")
  }
  
  #formChar <- paste("C",lis[1],"H",lis[2],"N",lis[3],"O",lis[4],"S",lis[5],"P",lis[6])
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







