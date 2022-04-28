getSAcontext <- function(SAcor, varType, altNuc, gen_cord, chromosome, 
                              strand, referenceDnaStringSet, location){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## Upstream proximal sequence
  upborder <- 55
  downborder <- 54
    
  ## Get reference sequence upstream of SA
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## Save ref sequence
  SA_sur <- sequence_range
  
  if(varType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }

  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                          indexCoordinate = SAcor,
                                          upRange = upborder, downRange = downborder, 
                                          strand = strand, referenceDnaStringSet)
  
  ## Get coordinate range
  donor_cords_up <- c((SAcor-upborder) : (SAcor+downborder))
  if(strand == "-1") donor_cords_up <- rev(c((SAcor-downborder) : (SAcor+upborder)))
  
  ## HZEI int is the same
  SA_surAlt <- SA_sur
  
  ## Calculate HZEI int in case variation lies within range
  if(gen_cord %in% donor_cords_up){
    
    SAHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_up %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) SAHBS_SNV[pos] <- altNuc
    if(varType == "DEL") SAHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") SAHBS_SNV[pos] <- paste0(SAHBS_SNV[pos], altNuc)
    SAHBS_SNV <-  paste(SAHBS_SNV, collapse = "")
    
    if(varType == "INS" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,1+insertion_length,
                                                   nchar(SAHBS_SNV))
    if(varType == "INS" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-insertion_length)
    if(varType == "DUP" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,2,
                                                   nchar(SAHBS_SNV))
    if(varType == "DUP" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-1)
    
    SA_surAlt <- SAHBS_SNV
    
  }
  
  return(c(SA_sur,SA_surAlt))
}
