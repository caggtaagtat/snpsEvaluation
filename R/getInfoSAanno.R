## Enter info about nearest SD
getInfoSAanno <- function(splicesites, varType, altNuc, gen_cord,
                          chromosome, strand, referenceDnaStringSet,
                          deletion_length, insertion_length){
  
  ## Convert strand annotation
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
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  listOfResults["SA MES ref"] <- as.numeric(calculateMaxEntScanScore(sequence_range, 3))
  
  ## Sequence range
  acceptor_cords <- c((SAcor-upborder) : (SAcor+downborder))
  if(strand == "-1") acceptor_cords <- rev(c((SAcor-downborder) : (SAcor+upborder)))
  
  varCoos <- gen_cord
  
  ## In case of deletion increase SA sequence by
  ## number of nt which overlapp with variation coordinates
  if(varType == "DEL"){
      
      ## check how many nucleotides overlap with SA sequence
      varCoos <- c(gen_cord:(gen_cord+deletion_length-1))
      varOverlapp <- sum(varCoos %in% acceptor_cords)
          
      if(location == "downstream")  downborder <- downborder+varOverlapp
      if(location == "upstream")  upborder <- upborder+varOverlapp
  }
  
  ## Get the surrounding reference sequence
  sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                         indexCoordinate = SAcor,
                                         upRange = upborder, downRange = downborder, 
                                         strand = strand, referenceDnaStringSet)
  
  ## If ref/alt SD has GT calculate MES
  listOfResults["SA MES delta"] <- 0
  listOfResults["SA MES alt"] <- listOfResults["SA MES ref"]
  
  ## In case genomic coordiante lies within 
  ## SA sequence re-calcualte MES score for alternative seq
  if(any(varCoos%in% acceptor_cords)){

    saMES_SNV <- strsplit(sequence_range, "")[[1]]
    
    ## Calcualte new MES
    pos <- which(acceptor_cords %in% gen_cord)
    if(varType == "SNV" | varType == "DUP" ) saMES_SNV[pos] <- altNuc
    if(varType == "DEL") saMES_SNV[pos:(pos+varOverlapp-1)] <- ""
    if(varType == "INS") saMES_SNV[pos] <- paste0(saMES_SNV[pos], altNuc)
    
    
    if(varType == "INS" & 
       location == "upstream") saMES_SNV <- saMES_SNV[1+insertion_length: nchar(saMES_SNV)]
    if(varType == "INS" & 
       location == "downstream") saMES_SNV <- saMES_SNV[1:(nchar(saMES_SNV)-insertion_length)]
    if(varType == "DUP" & 
       location == "upstream") saMES_SNV <- saMES_SNV[2: nchar(saMES_SNV)]
    if(varType == "DUP" & 
       location == "downstream") saMES_SNV <- saMES_SNV[1:(nchar(saMES_SNV)-1)]
    
    saMES_SNV <- paste(saMES_SNV, collapse = "")
    
    ## Save difference and alternative MES
    saMES_SNV <- as.numeric(calculateMaxEntScanScore(saMES_SNV, 3))
    if(length(saMES_SNV)==0) saMES_SNV <- -999
    saMES <- listOfResults$`SA MES ref`
    listOfResults["SA MES alt"] <- saMES_SNV
    listOfResults["SA MES delta"] <- round(saMES_SNV-saMES,digits=1)
    if(pos== 19 | pos== 20 | saMES_SNV == -999)  listOfResults["SA MES delta"] <- 9999
  }
  
  ## Save data fro output
  listOfResults$SAcor <- SAcor
  listOfResults$SAloc <- location
  
  return(listOfResults)
}
