location <- snvTable$SNV_localized_relative_to_SD[case]

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
  if(strand == "-1") donor_cords_up <- rev(donor_cords_up)
  
  #HZEI integral calculation
  SD_HZEI_up_alt <- SD_HZEI_up
  
  if(gen_cord %in% donor_cords_up){
    
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_up %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
    if(varType == "DEL") sdHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") sdHBS_SNV[pos] <- paste0(sdHBS_SNV[pos], altNuc)
    
    sdHBS_SNV <-  paste(sdHBS_SNV, collapse = "")
    SD_HZEI_up_alt <- calculateHZEIint(sdHBS_SNV)+0
    
  }
  
  
  #Get the respective genomic sequence downstream of the SD
  upborder <- 4
  downborder <- 63
  
  if(varType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }
  
  test2 <- Views(s[[as.character(chromosome)]], start=SDcor+upborder, end=SDcor+downborder)
  if(strand == -1)  test2 <- Views(s[[as.character(chromosome)]], start=SDcor-downborder, end=SDcor-upborder)
  
  test <- as.character(test2)
  sequence_range <-  as.character(test2)
  test3 <- data.frame(test2@ranges)
  
  #Minus strand sequence manipulation
  sequence_range[strand=="-1"] <- gsub("A","t",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("C","g",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("T","a",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("G","c",sequence_range[strand=="-1"])
  
  sequence_range[strand=="-1"] <- gsub("t","T",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("g","G",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("a","A",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- gsub("c","C",sequence_range[strand=="-1"])
  sequence_range[strand=="-1"] <- reverse(sequence_range[strand=="-1"])
  
  #Get coordinates
  donor_cords_down <- c((SDcor+4) : (SDcor+63))
  if(strand == "-1") donor_cords_down <- rev(c((SDcor-63) : (SDcor-4)))
  
  SD_seq_surround_down <- sequence_range
  
  #HZEI integral calculation
  SD_HZEI_down <- calculateHZEIint(sequence_range)+0
  SD_HZEI_down_alt <- calculateHZEIint(sequence_range)+0
  
  if(gen_cord %in% donor_cords_down){
    
    sdHBS_SNV <- strsplit(sequence_range, "")[[1]]
    
    pos <- which(donor_cords_down %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) sdHBS_SNV[pos] <- altNuc
    if(varType == "DEL") sdHBS_SNV[pos:(pos+deletion_length-1)] <- ""
    if(varType == "INS") sdHBS_SNV[pos] <- paste0(sdHBS_SNV[pos], altNuc)
    
    sdHBS_SNV[pos] <- altNuc
    sdHBS_SNV <-  paste(sdHBS_SNV, collapse = "")
    SD_HZEI_down_alt <- calculateHZEIint(sdHBS_SNV)+0
    
  }
  
  
  # Calculate SD SRE support
  res["SD SRE support"] <- SD_HZEI_up - SD_HZEI_down
  res["SD SRE support alternative"] <- SD_HZEI_up_alt - SD_HZEI_down_alt

  return(sdSREinfo)
}
