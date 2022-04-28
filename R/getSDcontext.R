getSDcontext <- function(SDcor, varType, altNuc, gen_cord, chromosome, 
                              strand, referenceDnaStringSet, location){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## Upstream proximal sequence
  upborder <- 57
  downborder <- 63
    
  ## Get reference sequence upstream of SD
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## SDve ref sequence
  SD_sur <- sequence_range
  
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
  
  ## HZEI int is the SDme
  SD_surAlt <- SD_sur
  
  ## Calculate HZEI int in case variation lies within range
  if(gen_cord %in% donor_cords_up){
    
    SDHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_up %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) SDHBS_SNV[pos] <- altNuc
    if(varType == "DEL") SDHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") SDHBS_SNV[pos] <- paste0(SDHBS_SNV[pos], altNuc)
    SDHBS_SNV <-  paste(SDHBS_SNV, collapse = "")
    
    if(varType == "INS" & 
       location == "upstream") SDHBS_SNV <- substr(SDHBS_SNV,1+insertion_length,
                                                   nchar(SDHBS_SNV))
    if(varType == "INS" & 
       location == "downstream") SDHBS_SNV <- substr(SDHBS_SNV,1,
                                                     nchar(SDHBS_SNV)-insertion_length)
    if(varType == "DUP" & 
       location == "upstream") SDHBS_SNV <- substr(SDHBS_SNV,2,
                                                   nchar(SDHBS_SNV))
    if(varType == "DUP" & 
       location == "downstream") SDHBS_SNV <- substr(SDHBS_SNV,1,
                                                     nchar(SDHBS_SNV)-1)
    
    SD_surAlt <- SDHBS_SNV
    
  }
  
  return(c(SD_sur,SD_surAlt))
}
