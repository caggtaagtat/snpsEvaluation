quickContext <- function(chromosome, gen_cord, strand, 
                         referenceDnaStringSet, varType, 
                         deletion_length){
  
    ## Get the surrounding reference sequence
    upborder <- 50
    downborder <- 50
    
    ## In case of an deletion, get longer sequence from reference
    ## Get the surrounding reference sequence
    sequence_range <- getReferenceSequence(chromosome = chromosome, 
                                           indexCoordinate = gen_cord,
                                           upRange = upborder, downRange = downborder, 
                                           strand = strand, referenceDnaStringSet)
    
    ## In case of a deletion widen the sequence range 
    ## in both directions by half of the deletion length
    if(varType == "DEL"){
        halfDelLength <- ceiling(deletion_length/2)
        downborder <- downborder + halfDelLength
        upborder <- upborder + halfDelLength
    }
    
    ## Get the surrounding reference sequence
    sequence_range2 <- getReferenceSequence(chromosome = chromosome, 
                                            indexCoordinate = gen_cord,
                                            upRange = upborder, downRange = downborder, 
                                            strand = strand, referenceDnaStringSet)
    
    ## Store sequence
    surroundingSeq  <- sequence_range
    alter <-  strsplit(sequence_range2, "")[[1]]
    
    ## Sequence range
    general_cords <- c((gen_cord-upborder) : (gen_cord+downborder))
    if(strand == "-1") general_cords <- rev(acceptor_cords)
    
    varPoint <- which(general_cords %in% gen_cord)
    
    ## Depending on variation type insert alternative sequence at position
    if(varType == "SNV" | varType == "DUP" ) alter[varPoint] <- altNuc
    if(varType == "DEL") alter[varPoint:(varPoint+deletion_length-1)] <- ""
    if(varType == "INS") alter[varPoint] <- paste0(alter[varPoint], altNuc)
    
    ## Save alternative sequence
    altSurroundingSeq <-  paste(alter, collapse="")
    
    ## Trim sequence in case of insertions or duplications
    if(varType == "INS"){
        if(insertion_length == 1){
            altSurroundingSeq <- substr(altSurroundingSeq, 2, nchar(altSurroundingSeq) )
        }else{
            his <- ceiling(insertion_length/2)
            altSurroundingSeq <- substr(altSurroundingSeq, 1+his, nchar(altSurroundingSeq)-his )
        }
    } 
    
    if(varType == "DUP") altSurroundingSeq <- substr(altSurroundingSeq, 2, nchar(altSurroundingSeq) )
    
    ## Calcualte Delta HZEI
    surroundingSeqHZEI <- calculateHZEIint(surroundingSeq)
    altSurroundingSeqHZEI <- calculateHZEIint(altSurroundingSeq)
    surroundingSeqHzeiDIFF <- altSurroundingSeqHZEI - surroundingSeqHZEI
    
  return(c(surroundingSeqHzeiDIFF, surroundingSeq, altSurroundingSeq ))
}
