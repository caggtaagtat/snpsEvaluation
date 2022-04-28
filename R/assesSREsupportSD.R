assesSREsupportSD <- function(SDcor, varType, altNuc, gen_cord, chromosome, 
                              strand, referenceDnaStringSet, location){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## Upstream proximal sequence
  ##
  ## Set borders 
  upborder <- 57
  downborder <- 2
    
  ## Get reference sequence upstream of SD
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## Get HZEI int of ref sequence
  SD_HZEI_up <- calculateHZEIint(sequence_range)+0

  if(varType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }

  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                          indexCoordinate = SDcor,
                                          upRange = upborder, downRange = downborder, 
                                          strand = strand, referenceDnaStringSet)
  
  ## Get coordinate range
  donor_cords_up <- c((SDcor-upborder) : (SDcor+downborder))
  if(strand == "-1") donor_cords_up <- rev(c((SDcor-downborder) : (SDcor+upborder)))
  
  ## HZEI int is the same
  SD_HZEI_up_alt <- SD_HZEI_up

  ## Calculate HZEI int in case variation lies within range
  if(gen_cord %in% donor_cords_up){
    
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_up %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
    if(varType == "DEL") sdHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") sdHBS_SNV[pos] <- paste0(sdHBS_SNV[pos], altNuc)
    sdHBS_SNV <-  paste(sdHBS_SNV, collapse = "")
    
    if(varType == "INS" & 
       location == "upstream") sdHBS_SNV <- substr(sdHBS_SNV,1+insertion_length,
                                                   nchar(sdHBS_SNV))
    if(varType == "INS" & 
       location == "downstream") sdHBS_SNV <- substr(sdHBS_SNV,1,
                                                     nchar(sdHBS_SNV)-insertion_length)
    if(varType == "DUP" & 
       location == "upstream") sdHBS_SNV <- substr(sdHBS_SNV,2,
                                                   nchar(sdHBS_SNV))
    if(varType == "DUP" & 
       location == "downstream") sdHBS_SNV <- substr(sdHBS_SNV,1,
                                                     nchar(sdHBS_SNV)-1)
    
    SD_HZEI_up_alt <- calculateHZEIint(sdHBS_SNV)+0

  }
  
  
  ## Get the respective genomic sequence 
  ## Downstream of the SD
  upborder <- -4
  downborder <- 63
  
  ## Get reference sequence upstream of SD
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## Get HZEI int of ref sequence
  SD_HZEI_down <- calculateHZEIint(sequence_range)+0

  if(varType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  ## Get coordinate range
  donor_cords_down <- c((SDcor-upborder) : (SDcor+downborder))
  if(strand == "-1") donor_cords_down <- rev(donor_cords_down)
  
  ## HZEI int is the same
  SD_HZEI_down_alt <- SD_HZEI_down
 
  ## In case variation lies within sequence range
  if(gen_cord %in% donor_cords_down){
    
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_down %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
    if(varType == "DEL") sdHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") sdHBS_SNV[pos] <- paste0(sdHBS_SNV[pos], altNuc)
    sdHBS_SNV <-  paste(sdHBS_SNV, collapse = "")
    
    if(varType == "INS" & 
       location == "upstream") sdHBS_SNV <- substr(sdHBS_SNV,1+insertion_length,
                                                   nchar(sdHBS_SNV))
    if(varType == "INS" & 
       location == "downstream") sdHBS_SNV <- substr(sdHBS_SNV,1,
                                                     nchar(sdHBS_SNV)-insertion_length)
    if(varType == "DUP" & 
       location == "upstream") sdHBS_SNV <- substr(sdHBS_SNV,2,
                                                   nchar(sdHBS_SNV))
    if(varType == "DUP" & 
       location == "downstream") sdHBS_SNV <- substr(sdHBS_SNV,1,
                                                     nchar(sdHBS_SNV)-1)
    
    ## Get alt HZEI integral
    sdHBS_SNV <- paste(sdHBS_SNV, collapse = "")
    SD_HZEI_down_alt <- calculateHZEIint(sdHBS_SNV)+0
  }
  
  
  # Calculate SD SRE support
  a <- SD_HZEI_up - SD_HZEI_down
  b <- SD_HZEI_up_alt - SD_HZEI_down_alt

  return(c(a,b))
}
