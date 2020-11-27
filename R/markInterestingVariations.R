markInterestingVariations <- function(SDdown, SDcor, SAdown, SAcor, sdHBS, saMES,
                                      HBSaltOver, HBSaltOverDiff, HBSaltOverDiffAno,
                                      MESaltOver, MESaltOverDiff, MESaltOverDiffAno,
                                      ){
  
  impactReport <- ""
  
  ###
  ## Changes for annotated SD
  ###
  
  ## Is GT dinucleotide of annotated SD affected
  if(SDdown == "yes") impactReport <- paste0("GT des annotierten SD betroffen (", SDcor, ")")

  ## Does Variation lead to new SD of over-average strength
  newSD <- HBSaltOver > 15 & HBSaltOverDiff > 2 & HBSaltOverDiffAno > -3.4
  if(newSD) impactReport <-  paste0("SD created/strengthened with/to ",HBSaltOver, " HBS. ",
                                    "Might compete with annotated SD (", SDcor,
                                " with ", sdHBS, " HBS)")
  
  ## Does Variation lead to new SA sites of over-average strength
  newSA <- MESaltOver > 8.4 & MESaltOverDiff > 1.8 & MESaltOverDiffAno > -3.3
  if(newSA) impactReport <-  paste0("SA created/strengthened with/to ",MESaltOver, " MaxEnt score. ",
                                    "Might compete with annotated SA (", SAcor,
                                    " with ",saMES, " MaxEnt score)")
  
 #MAxentAlt score treshold weniger strikt!
  
 ## SD SRE support changes drastically
 
  
  
  
  ###
  ## Changes for annotated SA
  ###
  
  ## Is AG dinucleotide of annotated SA affected
  if(SAdown == "yes") impactReport <- paste0("AG des annotierten SA betroffen (", SAcor, ")")
  
  
  ## Return dataframe at the end
  return(impactReport)
}
