markInterestingVariations <- function(SDdown, SDcor, SAdown, SAcor, sdHBS, saMES,
                                      HBSaltOver, HBSaltOverDiff, HBSaltOverDiffAno,
                                      MESaltOver, MESaltOverDiff, MESaltOverDiffAno,
                                      posToSD, posToSA, sdHBSalt, saMESalt, sdSREdiffperc,
                                      sdSRE, sdSREalt, saSRE, saSREalt, saSREdiffperc,
                                      dist_to_SD, dist_to_SA, GTcorSRE, GCcorSRE,
                                      AGcorSRE, details=TRUE){
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
  if(SDdiff < -2) impactReport <- paste0(impactReport,
                                             "HBS of annotated SD reduced from ",sdHBS,
                                         " to ", sdHBSalt, " (",
                                             SDcor, "). ")
  
  
  ## Is strength of annotated Acceptor affected
  SAdiff <- saMESalt - saMES
  if(SAdiff < -3.3) impactReport <- paste0(impactReport,
                                         "MaxEnt score of annotated SA reduced from ",saMES,
                                         " to ", saMESalt, " (",
                                         SAcor, "). ")
  
  
  
  
  ## Creats new donor
  if(posToSD == "upstream"){
    
    ## Does Variation lead to new SD upstream of the annotated one
    SDSREdiff <- GTcorSRE-sdSRE
    SDSREdiff <- SDSREdiff/sdSRE
    newSD <- HBSaltOverDiff > 2 & HBSaltOverDiffAno/sdHBS > -0.5 & SDSREdiff > -3 
    if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened with/to ",HBSaltOver, " HBS. ",
                                      "Might compete with ",SDposToSNV," annotated SD (", SDcor,
                                      " with ", sdHBS, " HBS). ")
  }
  
  if(posToSD == "downstream"){
    
    ## Does Variation lead to new SD upstream of the annotated one
    SDSREdiff <- GTcorSRE-sdSRE
    SDSREdiff <- SDSREdiff/sdSRE
    newSD <- HBSaltOverDiff > 2  & HBSaltOverDiffAno/sdHBS > -0.5 & SDSREdiff > -3 & dist_to_SD < 500
    if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened with/to ",HBSaltOver, " HBS. ",
                                      "Might compete with ",SDposToSNV," annotated SD (", SDcor,
                                      " with ", sdHBS, " HBS). ")
  }
  
  
  
  
  ## Creats new acceptor
  if(posToSA == "upstream"){
    SASREdiff <- AGcorSRE-saSRE
    SASREdiff <- SASREdiff/saSRE
    ## Does Variation lead to new SA site upstream of the annotated one
    newSA <- MESaltOverDiff > 3.3 & MESaltOverDiffAno/saMES > -0.5 & SASREdiff > -3
    if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened with/to ",MESaltOver, " MaxEnt score. ",
                                      "Might compete with ",SAposToSNV," annotated SA (", SAcor,
                                      " with ",saMES, " MaxEnt score). ")
   
  }
  
  if(posToSA == "downstream"){
    SASREdiff <- AGcorSRE-saSRE
    SASREdiff <- SASREdiff/saSRE
    ## Does Variation lead to new SA site downstream of the annotated one
    newSA <- MESaltOverDiff > 3.3 & MESaltOverDiffAno > -3.3 & SASREdiff > -3 & dist_to_SA < 500
    if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened with/to ",MESaltOver, " MaxEnt score. ",
                                      "Might compete with ",SAposToSNV," annotated SA (", SAcor,
                                      " with ",saMES, " MaxEnt score). ")
    
  }
  
  
  
  ## Looking at SRE support
  #sdSRE, sdSREalt, saSRE, saSREalt, sdSREdiffperc, saSREdiffperc
  
  sdSREdiff <- sdSREalt - sdSRE
  sdSREdiffperc <- abs(sdSREdiff/sdSRE*100)
  if(sdSREdiff < 0) sdSREdiffperc <- sdSREdiffperc*-1
  
  ## SD SRE support reduces by at least a fourth and makes the 
  ## alterative SRE support is below average
  sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE > 291.5 & sdSREdiff < -58
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SD (",
                                        SDcor,") decreases from over to below average. ")
  
  
  sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE < 291.5 & sdSREdiff < -58
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SD (",
                                        SDcor,") decreases. ")
  
  
  sdSREtest <- sdSREdiffperc < -25 &  sdHBSalt < 13.6 & sdSREdiff > -58
  
  if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SD (",
                                        SDcor,") decreases. ")


  
  
  
  ## SA SRE support reduces by at least a fourth and makes the 
  ## alterative SRE support is below average
  saSREdiff <- saSREalt - saSRE
  saSREdiffperc <- abs(saSREdiff/saSRE*100)
  if(saSREdiff < 0) saSREdiffperc <- saSREdiffperc*-1
  
  saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE > 377.2 & saSREdiff < -75
  
  if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SA (",
                                        SAcor,") decreases from over to below average. ")
  
  
  saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE < 377.2 & saSREdiff < -75
  
  if(saSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SA (",
                                        SAcor,") decreases. ")
  
  
  saSREtest <- saSREdiffperc < -25 &  saMESalt < 6.94 & saSREdiff < -75
  
  if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SA (",
                                        SAcor,") decreases. ")
  
  
  
  ## Return dataframe at the end
  if(details) return(impactReport)
  
  
  
  if(!details){
      
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
                                                 "GT of annotated SD affected")
      
      ## Is AG dinucleotide of annotated SA affected
      if(SAdown == "yes") impactReport <- paste0(impactReport,
                                                 "AG of annotated SA affected")
      
      ## Is strength of annotated Donor affected
      SDdiff <- sdHBSalt - sdHBS 
      if(SDdiff < -2) impactReport <- paste0(impactReport,
                                             "HBS of annotated SD reduced")
      
      
      ## Is strength of annotated Acceptor affected
      SAdiff <- saMESalt - saMES
      if(SAdiff < -3.3) impactReport <- paste0(impactReport,
                                               "MaxEnt score of annotated SA reduced")
      
      
      
      
      ## Creats new donor
      if(posToSD == "upstream"){
          
          ## Does Variation lead to new SD upstream of the annotated one
          SDSREdiff <- GTcorSRE-sdSRE
          SDSREdiff <- SDSREdiff/sdSRE
          newSD <- HBSaltOverDiff > 2 & HBSaltOverDiffAno/sdHBS > -0.5 & SDSREdiff > -3 
          if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened")
      }
      
      if(posToSD == "downstream"){
          
          ## Does Variation lead to new SD upstream of the annotated one
          SDSREdiff <- GTcorSRE-sdSRE
          SDSREdiff <- SDSREdiff/sdSRE
          newSD <- HBSaltOverDiff > 2  & HBSaltOverDiffAno/sdHBS > -0.5 & SDSREdiff > -3 & dist_to_SD < 500
          if(newSD) impactReport <-  paste0(impactReport,"SD created/strengthened")
      }
      
      
      
      
      ## Creats new acceptor
      if(posToSA == "upstream"){
          SASREdiff <- AGcorSRE-saSRE
          SASREdiff <- SASREdiff/saSRE
          ## Does Variation lead to new SA site upstream of the annotated one
          newSA <- MESaltOverDiff > 3.3 & MESaltOverDiffAno/saMES > -0.5 & SASREdiff > -3
          if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened")
          
      }
      
      if(posToSA == "downstream"){
          SASREdiff <- AGcorSRE-saSRE
          SASREdiff <- SASREdiff/saSRE
          ## Does Variation lead to new SA site downstream of the annotated one
          newSA <- MESaltOverDiff > 3.3 & MESaltOverDiffAno > -3.3 & SASREdiff > -3 & dist_to_SA < 500
          if(newSA) impactReport <-  paste0(impactReport,"SA created/strengthened")
          
      }
      
      
      
      ## Looking at SRE support
      #sdSRE, sdSREalt, saSRE, saSREalt, sdSREdiffperc, saSREdiffperc
      
      sdSREdiff <- sdSREalt - sdSRE
      sdSREdiffperc <- abs(sdSREdiff/sdSRE*100)
      if(sdSREdiff < 0) sdSREdiffperc <- sdSREdiffperc*-1
      
      ## SD SRE support reduces by at least a fourth and makes the 
      ## alterative SRE support is below average
      sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE > 291.5 & sdSREdiff < -58
      
      if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SD",
                                            " decreases from over to below average")
      
      
      sdSREtest <- sdSREdiffperc < -25 &  sdSREalt < 291.5 &  sdSRE < 291.5 & sdSREdiff < -58
      
      if(sdSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SD decreases")
      
      
      sdSREtest <- sdSREdiffperc < -25 &  sdHBSalt < 13.6 & sdSREdiff > -58
      
      if(sdSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SD",
                                            " decreases")
      
      
      
      
      
      ## SA SRE support reduces by at least a fourth and makes the 
      ## alterative SRE support is below average
      saSREdiff <- saSREalt - saSRE
      saSREdiffperc <- abs(saSREdiff/saSRE*100)
      if(saSREdiff < 0) saSREdiffperc <- saSREdiffperc*-1
      
      saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE > 377.2 & saSREdiff < -75
      
      if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of annotated SA",
                                            " decreases from over to below average")
      
      
      saSREtest <- saSREdiffperc < -25 &  saSREalt < 377.2 &  saSRE < 377.2 & saSREdiff < -75
      
      if(saSREtest) impactReport <-  paste0(impactReport,"Below average SRE support of annotated SA",
                                            " decreases")
      
      
      saSREtest <- saSREdiffperc < -25 &  saMESalt < 6.94 & saSREdiff < -75
      
      if(saSREtest) impactReport <-  paste0(impactReport,"SRE support of weak annotated SA ",
                                            " decreases")
      
      
      
      ## Return dataframe at the end
      return(impactReport)
  }
}
