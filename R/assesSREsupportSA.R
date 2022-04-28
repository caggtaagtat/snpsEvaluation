assesSREsupportSA <- function(SAcor, variationType, altNuc, gen_cord, chromosome, 
                              strand, referenceDnaStringSet, location){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  ## Upstream proximal sequence
  ##
  ## Set borders 
  upborder <- 55
  downborder <- 4
    
  ## Get reference sequence upstream of SA
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## Get HZEI int of ref sequence
  SA_HZEI_up <- calculateHZEIint(sequence_range)+0
  
  if(variationType == "DEL"){
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
  SA_HZEI_up_alt <- SA_HZEI_up
  
  ## Calculate HZEI int in case variation lies within range
  if(gen_cord %in% donor_cords_up){
    
    SAHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_up %in% gen_cord)
    if(variationType == "SNV" | variationType == "DUP" ) SAHBS_SNV[pos] <- altNuc
    if(variationType == "DEL") SAHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(variationType == "INS") SAHBS_SNV[pos] <- paste0(SAHBS_SNV[pos], altNuc)
    SAHBS_SNV <-  paste(SAHBS_SNV, collapse = "")
    
    if(variationType == "INS" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,1+insertion_length,
                                                   nchar(SAHBS_SNV))
    if(variationType == "INS" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-insertion_length)
    if(variationType == "DUP" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,2,
                                                   nchar(SAHBS_SNV))
    if(variationType == "DUP" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-1)
    
    SA_HZEI_up_alt <- calculateHZEIint(SAHBS_SNV)+0
    
  }
  
  
  ## Get the respective genomic sequence 
  ## Downstream of the SA
  upborder <- 5
  downborder <- 54
  
  ## Get reference sequence upstream of SA
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  ## Get HZEI int of ref sequence
  SA_HZEI_down <- calculateHZEIint(sequence_range)+0

  if(variationType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  ## Get coordinate range
  donor_cords_down <- c((SAcor-upborder) : (SAcor+downborder))
  if(strand == "-1") donor_cords_down <- rev(donor_cords_down)
  
  ## HZEI int is the same
  SA_HZEI_down_alt <- SA_HZEI_down
 
  ## In case variation lies within sequence range
  if(gen_cord %in% donor_cords_down){
    
    SAHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_down %in% gen_cord)
    if(variationType == "SNV" | variationType == "DUP" ) SAHBS_SNV[pos] <- altNuc
    if(variationType == "DEL") SAHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(variationType == "INS") SAHBS_SNV[pos] <- paste0(SAHBS_SNV[pos], altNuc)
    SAHBS_SNV <-  paste(SAHBS_SNV, collapse = "")
    
    if(variationType == "INS" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,1+insertion_length,
                                                   nchar(SAHBS_SNV))
    if(variationType == "INS" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-insertion_length)
    if(variationType == "DUP" & 
       location == "upstream") SAHBS_SNV <- substr(SAHBS_SNV,2,
                                                   nchar(SAHBS_SNV))
    if(variationType == "DUP" & 
       location == "downstream") SAHBS_SNV <- substr(SAHBS_SNV,1,
                                                     nchar(SAHBS_SNV)-1)
    
    ## Get alt HZEI integral
    SAHBS_SNV <- paste(SAHBS_SNV, collapse = "")
    SA_HZEI_down_alt <- calculateHZEIint(SAHBS_SNV)+0
    
  }
  
  
  # Calculate SA SRE support
  a <-  SA_HZEI_down - SA_HZEI_up
  b <-  SA_HZEI_down_alt - SA_HZEI_up_alt

  return(c(a,b))
}
