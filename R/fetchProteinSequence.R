fetchProteinSequence <- function( uniprotSpeciesName = "Homo sapiens",  proteinAccession = "P02144"){
  
  tmp <- availableUniprotSpecies()
  
  if (any(tmp$`Species name` == uniprotSpeciesName)) taxId <- as.integer(tmp$`taxon ID`[tmp$`Species name` == uniprotSpeciesName])
  #else generate warning and quit
  
  tmp <- UniProt.ws(taxId)
  
  proteinSequence <- select(tmp, proteinAccession, "SEQUENCE", "UNIPROTKB")[,2]
  
  return(proteinSequence)
}