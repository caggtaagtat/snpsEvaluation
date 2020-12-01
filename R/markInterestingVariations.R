markInterestingVariations <- function(SDdown, SDcor, SAdown, SAcor, sdHBS, saMES,
                                      HBSaltOver, HBSaltOverDiff, HBSaltOverDiffAno,
                                      MESaltOver, MESaltOverDiff, MESaltOverDiffAno,
                                      posToSD, posToSA, sdHBSalt, saMESalt, sdSREdiffperc,
                                      sdSRE, sdSREalt, saSRE, saSREalt, saSREdiffperc){
  
  ##Invert location
  if(posToSD == "upstream") SDposToSNV <- "downstream"
  if(posToSD == "downstream") SDposToSNV <- "upstream"
  
  if(posToSA == "upstream") SAposToSNV <- "downstream"
  if(posToSA == "downstream") SAposToSNV <- "upstream"
  
  impactReport <- ""
  
  ###
  ## Changes for annotated SD
  ###
  
  ## Is GT dinucleotide of annotated SD affected
  if(SDdown == "yes") impactReport <- paste0(impactReport,
                                             "GT of annotated SD affected (",
                                             SDcor, "). ")

  ## Is AG dinucleotide of annotated SA affected
  if(SAdown == "yes") impactReport <- paste0(impactReport,
                                             "AG of annotated SA affected (",
                                             SAcor, "). ")
  
  ## Is strength of annotated Donor affected
  SDdiff <- sdHBSalt - sdHBS
  if(SDdiff < -3.4) impactReport <- paste0(impactReport,
                                             "HBS of annotated SD reduced by ",SDdiff,
                                         " to ", sdHBSalt, " (",
                                             SDcor, "). ")
  
  
  ## Is strength of annotated Acceptor affected
  SAdiff <- saMESalt - saMES
  if(SAdiff < -3.4) impactReport <- paste0(impactReport,
                                         "MaxEnt score of annotated SA reduced by ",SAdiff,
                                         " to ", saMESalt, " (",
                                         SAcor, "). ")
  
  
  
  
  ## Creats new donor
  if(posToSD == "upstream"){
    
    ## Does Variation lead to new SD upstream of the annotated one
    newSD <- HBSaltOverDiff > 3.4 & HBSaltOverDiffAno > -3.4 & HBSaltOverDiffAno != 0
    if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened with/to ",HBSaltOver, " HBS. ",
                                      "Might compete with ",SDposToSNV," annotated SD (", SDcor,
                                      " with ", sdHBS, " HBS). ")
  }
  
  if(posToSD == "downstream"){
    
    ## Does Variation lead to new SD upstream of the annotated one
    newSD <- HBSaltOverDiff > 3.4 & HBSaltOverDiffAno > -1 & HBSaltOverDiffAno != 0
    if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened with/to ",HBSaltOver, " HBS. ",
                                      "Might compete with ",SDposToSNV," annotated SD (", SDcor,
                                      " with ", sdHBS, " HBS). ")
  }
  
  
  
  
  
  
  ## Creats new acceptor
  if(posToSA == "upstream"){
    
    ## Does Variation lead to new SA site upstream of the annotated one
    newSA <- MESaltOverDiff > 3.4 & MESaltOverDiffAno > -1 & MESaltOverDiffAno != 0
    if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened with/to ",MESaltOver, " MaxEnt score. ",
                                      "Might compete with ",SAposToSNV," annotated SA (", SAcor,
                                      " with ",saMES, " MaxEnt score). ")
   
  }
  
  if(posToSA == "downstream"){
    
    ## Does Variation lead to new SA site downstream of the annotated one
    newSA <- MESaltOverDiff > 3.4 & MESaltOverDiffAno > -3.4 & MESaltOverDiffAno != 0
    if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened with/to ",MESaltOver, " MaxEnt score. ",
                                      "Might compete with ",SAposToSNV," annotated SA (", SAcor,
                                      " with ",saMES, " MaxEnt score). ")
    
  }
  
  
  
  ## Looking at SRE support
  #sdSRE, sdSREalt, saSRE, saSREalt, sdSREdiffperc, saSREdiffperc
  
  ## SD SRE support reduces by at least a fourth and makes the 
  ## alterative SRE support is below average
  sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE > 291.5
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SD (",
                                        SDcor,") decreases from over to below average. ")
  
  
  sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE < 291.5
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SD (",
                                        SDcor,") decreases. ")
  
  
  sdSREtest <- sdSREdiffperc < -25 &  sdHBSalt < 13.6
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SD (",
                                        SDcor,") decreases. ")


  
  
  
  ## SA SRE support reduces by at least a fourth and makes the 
  ## alterative SRE support is below average
  saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE > 377.2
  
  if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SA (",
                                        SAcor,") decreases from over to below average. ")
  
  
  saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE < 377.2
  
  if(saSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SA (",
                                        SAcor,") decreases. ")
  
  
  saSREtest <- saSREdiffperc < -25 &  saMESalt < 6.94
  
  if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SA (",
                                        SAcor,") decreases. ")
  
  
  
  ## Return dataframe at the end
  return(impactReport)
}
