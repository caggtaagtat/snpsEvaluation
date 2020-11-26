## Define input variables for manual testing
# variationType <- varType
# alt_nuc <- altNuc
# genomic_coordinate <- gen_cord

## Enter info about nearest SD
getInfoSAanno <- function(splicesites, variationType, alt_nuc, genomic_coordinate,
                          chromosome, strand, referenceDnaStringSet){
  
  #Convert strand annotation
  if(strand== "+") strand <- "1"
  if(strand== "-") strand <- "-1"
  
  listOfResults <- list()
  
  ## Location of variation
  location <- "downstream"
  if(splicesites$distance[splicesites$type == "SA"][1] < 0) location <- "upstream"
  
  listOfResults["SNP distance to SA NR"] <- splicesites$distance_abs[splicesites$type == "SA"][1]
  SAcor <- splicesites$coordinate[splicesites$type == "SA"][1]
  
  ## Get Sequence of annotated Acceptor with and without the variation
  upborder <- 20
  downborder <- 2
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  listOfResults["SA MES ref"] <- as.numeric(calculateMaxEntScanScore(sequence_range, 3))
  
  if(variationType == "DEL"){
    if(location == "downstream")  downborder <- downborder+deletion_length
    if(location == "upstream")  upborder <- upborder+deletion_length
  }
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SDcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  ## Sequence range
  acceptor_cords <- c((SAcor-upborder) : (SAcor+downborder))
  if(strand == "-1") acceptor_cords <- rev(acceptor_cords)
  
  ## If ref/alt SD has GT calculate MES
  listOfResults["SA MES delta"] <- 0
  listOfResults["SA MES alt"] <- listOfResults["SA MES ref"]
  
  ## In case genomic coordiante lies within 
  ## SA sequence re-calcualte MES score for alternative seq
  if(genomic_coordinate %in% acceptor_cords){

    saMES_SNV <- strsplit(sequence_range, "")[[1]]
    
    ## Calcualte new MES
    pos <- which(acceptor_cords %in% genomic_coordinate)
    if(variationType == "SNV" | variationType == "DUP" ) saMES_SNV[pos] <- alt_nuc
    if(variationType == "DEL") saMES_SNV[pos:(pos+deletion_length-1)] <- ""
    if(variationType == "INS") saMES_SNV[pos] <- paste0(saMES_SNV[pos], alt_nuc)
    saMES_SNV <- paste(saMES_SNV, collapse = "")
    
    
    if(variationType == "INS" & 
       location == "upstream") saMES_SNV <- substr(saMES_SNV,1+insertion_length,
                                                   nchar(saMES_SNV))
    if(variationType == "INS" & 
       location == "downstream") saMES_SNV <- substr(saMES_SNV,1,
                                                     nchar(saMES_SNV)-insertion_length)
    if(variationType == "DUP" & 
       location == "upstream") saMES_SNV <- substr(saMES_SNV,2,
                                                   nchar(saMES_SNV))
    if(variationType == "DUP" & 
       location == "downstream") saMES_SNV <- substr(saMES_SNV,1,
                                                     nchar(saMES_SNV)-1)
    
    ## Save difference and alternative MES
    saMES_SNV <- as.numeric(calculateMaxEntScanScore(saMES_SNV, 3))
    if(length(saMES_SNV)==0) saMES_SNV <- -999
    saMES <- listOfResults["SA MES ref"] 
    listOfResults["SA MES alt"] <- saMES_SNV
    listOfResults["SA MES delta"] <- round(saMES_SNV-saMES,digits=1)
    if(pos== 19 | pos== 20 | saMES_SNV == -999)  listOfResults["SA MES delta"] <- 9999
  }
  
  listOfResults$SAcor <- SAcor
  listOfResults$SAloc <- location
  return(listOfResults)
}
