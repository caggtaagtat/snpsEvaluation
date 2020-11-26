genomicPosAndNTs <- function(anno, referenceDnaStringSet, transcriptID, transCur){
  
  deletion_length <- 0
  insertion_length <- 0
  
  #Get variation
  if(str_detect(anno, "dup")==F){
    find_arrow <- regexpr(">", anno)[[1]]
    from_coding <- substr(anno,3, find_arrow-2 )
    refNuc <- substr(anno,find_arrow-1,find_arrow-1 )
    altNuc <- substr(anno,find_arrow+1,find_arrow+1 )
    if(str_detect(anno, "del")==F){
      
      if(str_detect(anno, "ins")==T){
        altNuc <- substr(anno, nchar(anno), nchar(anno))
        refNuc <- ""
        find_arrow <- regexpr("_", anno)[[1]]
        if(find_arrow != -1){
          find_del <- regexpr("ins", anno)[[1]]
          first_coordinate <- substr(anno,3, find_arrow-1 )
          second_coordinate <- substr(anno, find_arrow+1, find_del-1)
          insertion_length <- length(c(first_coordinate : second_coordinate))
          pseudoVariation <- paste0("c.",first_coordinate,"A>T")
        }else{
        pseudoVariation <- paste0(strsplit(anno, "ins")[[1]][[1]], "A>T")
        insertion_length <- 1
        }
        
      }else{
        pseudoVariation <- anno
      }
      
      results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID,
                                         pseudoVariation , ntWindow=20, transCur,
                                         gene2transcript=gene2transcript)
      
      replace_coord <- results$genomicCoordinate
      
    }
    
    if(str_detect(anno, "del")==T){
      
      refNuc <- substr(anno, nchar(anno), nchar(anno))
      altNuc <- ""
      find_arrow <- regexpr("_", anno)[[1]]
      if(find_arrow != -1){
        find_del <- regexpr("del", anno)[[1]]
        first_coordinate <- substr(anno,3, find_arrow-1 )
        second_coordinate <- substr(anno, find_arrow+1, find_del-1)
        deletion_length <- length(c(first_coordinate : second_coordinate))
        pseudoVariation <- paste0("c.",first_coordinate,"A>T")
      }else{
        pseudoVariation <- paste0(strsplit(anno, "del")[[1]][[1]], "A>T")
        deletion_length <- 1
        
      }
      
      
      results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID,
                                         pseudoVariation, ntWindow=20, transCur,
                                         gene2transcript=gene2transcript)
      
      replace_coord <- results$genomicCoordinate
    }
    
  }else{
    refNuc <- substr(anno, nchar(anno), nchar(anno))
    altNuc <- paste0(refNuc, refNuc)
    pseudoVariation <- paste0(strsplit(anno, "dup")[[1]][[1]], "A>T")
    
    results <- getSeqInfoFromVariation(referenceDnaStringSet, transcriptID,
                                       pseudoVariation, ntWindow=20, transCur,
                                       gene2transcript=gene2transcript)
    
    replace_coord <- results$genomicCoordinate
  }
  
  return(paste(replace_coord, refNuc, altNuc, deletion_length, insertion_length))
}
