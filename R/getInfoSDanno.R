## Enter info about nearest SD
getInfoSDanno <- function(splicesites, varType, altNuc, gen_cord,
                          chromosome, strand, referenceDnaStringSet,deletion_length,
                          insertion_length){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  listOfResults <- list()
  
  ## Location of variation
  location <- "downstream"
  if(splicesites$distance[splicesites$type == "SD"][1] < 0) location <- "upstream"
  
  listOfResults["SNP distance to SD NR"] <- splicesites$distance_abs[splicesites$type == "SD"][1]
  SDcor <- splicesites$coordinate[splicesites$type == "SD"][1]
  
  ## Get Sequence of annotated Donor with and without the variation
  upborder <- 2
  downborder <- 8
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                          indexCoordinate = SDcor,
                                          upRange = upborder, downRange = downborder, 
                                          strand = strand, referenceDnaStringSet)
  
  if(substr(sequence_range, 4,5) == "GT") listOfResults["SD HBS ref"] <- hbg$hbs[hbg$seq == sequence_range] 
  if(substr(sequence_range, 4,5) != "GT") listOfResults["SD HBS ref"] <- 999
  
  ## Sequence range
  donor_cords <- c((SDcor-upborder) : (SDcor+downborder))
  if(strand == "-1") donor_cords <- rev(c((SDcor-downborder) : (SDcor+upborder)))
  
  ## In case of deletion increase SD sequence by
  ## number of nt which overlapp with variation coordinates
  if(varType == "DEL"){
      ## check how many nucleotides overlap with SA sequence
      varCoos <- c(gen_cord:(gen_cord+deletion_length-1))
      varOverlapp <- sum(varCoos %in% donor_cords)
      
      if(location == "downstream")  downborder <- downborder+deletion_length
      if(location == "upstream")  upborder <- upborder+deletion_length
  }

  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)

  ## If ref/alt SD has GT calculate HBS
  listOfResults["SD HBS delta"] <- 0
  listOfResults["SD HBS alt"] <- listOfResults["SD HBS ref"]
  
  ## In case genomic coordiante lies within 
  ## SD sequence re-calcualte HBS score for alternative seq
  if(gen_cord %in% donor_cords){
    
    sdHBS <- listOfResults$'SD HBS ref'
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    ## Calcualte new HBS
    pos <- which(donor_cords %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
    if(varType == "DEL") sdHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") sdHBS_SNV[pos] <- paste0(sdHBS_SNV[pos], altNuc)
    sdHBS_SNV <- paste(sdHBS_SNV, collapse = "")
    
    if(varType == "INS" & 
       location == "upstream") sdHBS_SNV <- sdHBS_SNV[1+insertion_length: nchar(sdHBS_SNV)]
    if(varType == "INS" & 
       location == "downstream") sdHBS_SNV <- sdHBS_SNV[1:(nchar(sdHBS_SNV)-insertion_length)]
    if(varType == "DUP" & 
       location == "upstream") sdHBS_SNV <- sdHBS_SNV[2: nchar(sdHBS_SNV)]
    if(varType == "DUP" & 
       location == "downstream") sdHBS_SNV <- sdHBS_SNV[1:(nchar(sdHBS_SNV)-1)]
    
    sdHBS_SNV <- hbg$hbs[hbg$seq == sdHBS_SNV] 
    
    ## Save difference and alternative HBS
    if(length(sdHBS_SNV)==0) sdHBS_SNV <- -999
    listOfResults["SD HBS alt"] <- sdHBS_SNV
    listOfResults["SD HBS delta"] <- round(sdHBS_SNV-sdHBS,digits=1)
    if(pos==4 | pos==5 | sdHBS_SNV == -999 | sdHBS == 999 )  listOfResults["SD HBS delta"] <- 9999
  }
  
  listOfResults$SDcor <- SDcor
  listOfResults$SDloc <- location
  return(listOfResults)
}
