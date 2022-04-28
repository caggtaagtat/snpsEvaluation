genomicPosAndNTs <- function(anno, referenceDnaStringSet, transcript, transCur){
  
  ## Default output value
  deletion_length <- 0
  insertion_length <- 0
  
  ## In case variation is not a duplication
  if(str_detect(anno, "dup")==F){
    
    ## Get postion of > character and nucleotides
    find_arrow <- regexpr(">", anno)[[1]]
    from_coding <- substr(anno,3, find_arrow-2 )
    refNuc <- substr(anno,find_arrow-1,find_arrow-1 )
    altNuc <- substr(anno,find_arrow+1,find_arrow+1 )
    
    ## In case variaton is not deletion
    if(str_detect(anno, "del")==F){
      
      ## If variations is an insertion
      if(str_detect(anno, "ins")==T){
        
        ## if there is only a single nucleotide insertion
        ## Get inserted nucleotide as alternative
        altNuc <- substr(anno, nchar(anno), nchar(anno))
        refNuc <- ""
        find_arrow <- regexpr("_", anno)[[1]]
        
        ## In case the insertion annotation has a "_" character
        if(find_arrow != -1){
          
          ## Find "ins" character position and insertion coordinates
          find_ins <- regexpr("ins", anno)[[1]]
          first_coordinate <- substr(anno,3, find_arrow-1 )
          second_coordinate <- substr(anno, find_arrow+1, find_ins-1)
          
          ## Get insertion length and insertion sequence
          insertion_length <- length(c(first_coordinate : second_coordinate))
          altNuc <- substr(anno, find_ins+3, nchar(anno))
          
          ## Create pseudo annotation for coordiante calculations
          pseudoVariation <- paste0("c.",first_coordinate,"A>T")
        
        ## In case insertion annotation has NOT a "_" character  
        }else{
          
        ## Create pseudo annotation for coordiante calculations
        pseudoVariation <- paste0(strsplit(anno, "ins")[[1]][[1]], "A>T")
        insertion_length <- 1
        }
      
      ## If variations is NOT an insertion and NOT a deletion  
      }else{
        
        ## Use real annotation in pseudovariable
        pseudoVariation <- anno
      }
      
      ## Convert coding position into genomic position
      results <- getSeqInfoFromVariation(referenceDnaStringSet, transcript,
                                         pseudoVariation , ntWindow=20, transCur,
                                         gene2transcript=gene2transcript)
      ## Save retrieve coordinate
      replace_coord <- results$genomicCoordinate
      
    }
    
    ## In case variaton is a deletion
    if(str_detect(anno, "del")==T){
      
      ## If its a single nucleotide deletion get the reference nt
      refNuc <- substr(anno, nchar(anno), nchar(anno))
      altNuc <- ""
      find_arrow <- regexpr("_", anno)[[1]]
      
      ## In case deletion is longer than 1 NT
      if(find_arrow != -1){
        
        ## Find "del" character in annotation and the coordiantes and length
        find_del <- regexpr("del", anno)[[1]]
        first_coordinate <- substr(anno,3, find_arrow-1 )
        second_coordinate <- substr(anno, find_arrow+1, find_del-1)
        deletion_length <- length(c(first_coordinate : second_coordinate))
        
        ## Use first coodinate in pseuovariation annotation
        pseudoVariation <- paste0("c.",first_coordinate,"A>T")
        
      ## In case deletion is NOT longer than 1 NT
      }else{
        
        ## Use stated coordinate as pseudoannotation
        pseudoVariation <- paste0(strsplit(anno, "del")[[1]][[1]], "A>T")
        deletion_length <- 1
      }
      
      ## Convert translational position into genomic
      results <- getSeqInfoFromVariation(referenceDnaStringSet, transcript,
                                         pseudoVariation, ntWindow=20, transCur,
                                         gene2transcript=gene2transcript)
      ## Save results
      replace_coord <- results$genomicCoordinate
    }
    
  ## If variation is a duplication  
  }else{
    
    ## Get alt and ref NT from annotation
    refNuc <- substr(anno, nchar(anno), nchar(anno))
    altNuc <- paste0(refNuc, refNuc)
    
    ## Create pseudo-annotation for coordiante conversion
    pseudoVariation <- paste0(strsplit(anno, "dup")[[1]][[1]], "A>T")
    
    ## Convert translational position into genomic
    results <- getSeqInfoFromVariation(referenceDnaStringSet, transcript,
                                       pseudoVariation, ntWindow=20, transCur,
                                       gene2transcript=gene2transcript)
    ## Save results
    replace_coord <- results$genomicCoordinate
  }
  
  ## Return all saved parameters
  return(paste(replace_coord, refNuc, altNuc, deletion_length, insertion_length))
}
