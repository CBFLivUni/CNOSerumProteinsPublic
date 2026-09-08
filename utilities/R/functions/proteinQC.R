

GetSerumQCProteins <- function(serumQCData) {
  
  # Get serum proteins
  serumProtsCheck <- apply(serumQCData, 1, function(prot) all(!is.na(prot) & !is.nan(prot)))
  serumProts <- rownames(serumQCData)[serumProtsCheck]; length(serumProts)
  
  # Return
  return(serumProts)
  
}
