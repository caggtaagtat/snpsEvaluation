createOutputDocu <- function(){
  
  columns <- c("gene",
               "c.DNA",
               "Delta HZEI integral",
               "Delta MaxEnt anno SA",
               "annoSA AG dinucleotide affected",
               "SD HBS diff",
               "annoSD GT dinucleotide affected",
               "SNP distance to SA nt",
               "SNP distance to SD nt",
               "transcript context",
               "SA SRE support",
               "SA SRE support alternative",
               "Delta SA SRE support",
               "Delta SA SRE support percent",
               "SD SRE support",
               "SD SRE support alternative",
               "Delta SD SRE support",
               "Delta SD SRE support percent",
               "HBS alt overlapping SDs",
               "HBS diff overlapping SDs",
               "HBS diff overlapping SDs vs annotated",
               "MaxEnt alt overlapping SAs",
               "MaxEnt diff overlapping SAs",
               "MaxEnt diff overlapping SA vs annotated",
               "Maxent alt overlapping SD GC sites",
               "Maxent diff overlapping SD GC sites",
               "SD HBS ref",
               "SD HBS alt",
               "SA MaxEnt ref",
               "SA MaxEnt alt",
               "SNV localized relative to SD",
               "SNV localized relative to SA",
               "length of next exon",
               "variation type",
               "second up ss strength",
               "second down ss strength",
               "gene strand",
               "chromosome",
               "ENSEMBL transcriptID",
               "genome position",
               "alt",
               "Exon sequence",
               "SD sequence surrounding",
               "SD sequence surrounding with SNP",
               "SA sequence surrounding",
               "SA sequence surrounding with SNP",
               "Potentially interesting",
               "ID")
  
  documentation <- c(
       "gene name",  
       
       "Annotation der Variation",
       
       "Differenz im HZEI-integral durch die variation",
       
       "Änderungen im MaxEnt score des annotierten Akzeptor",
       
       "Ist das AG Dinucleotid des annotierten Akzeptor betroffen?",
                                                                                                                                                                                                                                                                                                               
       paste0("Änderungen im HBond score des annotierten Donors. Ist 9999,", 
              "falls das GT des Donors durch die Variation gestört wird"),                                                                                                                                                       
       
       "Ist das GT Dinucleotid des annotierten Donors betroffen?",  
       
       "Entfernung der variation zum nächsten annotierten SD", 
       
       "Entfernung der Variation zum nächsten annotierten SA", 
       
       "Ist die Variation exonisch oder intronisch im Transkript",
       
       paste0("HZEI-Integral downstream minus HZEI-Integral upstream des", 
       "nächstgelegenen SA in Reference-Sequzenz"),   
       
       paste0("HZEI-Integral downstream minus HZEI-Integral upstream des", 
              "nächstgelegenen SA in Reference-Sequzenz mit der Variation"),                                                                                                                                                                                                                     
       
       paste0("Differenz im HZEI-Integral Gewicht des SA mit und ohne",  
       "Variation"),    
       
       paste0("(Differenz im HZEI-Integral Gewicht des SA) geteilt durch", 
              "SA SRE support in %"),              
       
       paste0("HZEI-Integral upstream minus HZEI-Integral downstream des",  
       "nächstgelegenen SD in Reference-Sequzenz"),          
       
       paste0("HZEI-Integral upstream minus HZEI-Integral downstream des", 
              "nächstgelegenen SD mit der Variation"),     
       
       "Differenz im HZEI-Integral Gewicht des SD", 
       
       paste0("(Differenz im HZEI-Integral Gewicht des SD) geteilt durch", 
       "SD SRE support in %"),                         
       
       paste0("Maximaler HBS von überlappenden SD-Sequenzen in",  
              "Referenz-Sequenz mit eingebauter Variation"),    
       
       paste0("Differenz des maximalen HBS von überlappenden",  
       "SD-Sequenzen. Unterschied nach Einbauen der Variation",
              "in Referenz-Sequenz"),     
       
       paste0("Differenz des maximalen HBS von überlappenden",  
       "SD-Sequenzen zu dem HBS des nächstegelegenen Donors.",  
       "Unterschied nach Einbauen der Variation
       in Referenz-Sequenz"), 
       
       paste0("Maximaler MaxEnt score von überlappenden SA-Sequenzen", 
       "in Referenz-Sequenz mit eingebauter Variation"),     
       
       paste0("Differenz des maximalen MaxEnt scores  von überlappenden", 
              "SA-Sequenzen. Unterschied nach Einbauen der Variation", 
               "in Referenz-Sequenz"), 
       
       paste0("Differenz des maximalen MaxEnt scores von überlappenden",  
              "SA-Sequenzen zu dem MaxEnt scores des nächstegelegenen",  
       "Akzeptors Unterschied nach Einbauen der Variation", 
              "in Referenz-Sequenz"), 
       
       paste0("Maximaler MaxEnt score von überlappenden SD-Sequenzen", 
       "mit einem GC-Dinukleotid auf SD Position +/+ in ", 
              "Referenz-Sequenz mit eingebauter Variation"),    
       
       paste0("Differenz des maximalen MaxEnt score  von überlappenden", 
       "SD-Sequenzen mit einem GC-Dinukleotid auf SD Position +/+.", 
              "Unterschied nach Einbauen der Variation in Referenz-Sequenz"),  
       
       "HBS des annotierten SD in Referenzgenom",       
       
       paste0("HBS des annotierten SD mit der entsprechenden Variation", 
       "(-999 falls resultierende SD mit der Variation kein GT mehr hat)"), 
       
       "MaxEnt score des Akzeptors im Referenzgenom",  
       
       "MaxEnt score des Akzeptors mit der enstprechenden Variation", 
       
       paste0("Liegt der Variation upstream oder downstream zum nächsten", 
       "annotierten SD"),             
       
       paste0("Liegt der Variation upstream oder downstream zum nächsten", 
              "annotierten SA"),                                         
       
       "Länge des nächst-gelegenen Exons",                   
       
       "Typ der Variation (Variation, DEL, DUP oder INS)", 
       
       paste0("Splice site Stärke der übernächsten splice site upstream der", 
       " Variation"),                                           
       
       paste0("Splice site Stärke der übernächsten splice site downstream der", 
              "Variation"),                                       
       
       "gene_strand",   
       
       "chromosome", 
       
       "ENSEMBL_transcriptID",
       
       "genome_position",   
       
       "alt",            
       
       "Exon_sequence" , 
       
       paste0("Umgebendes Sequenz-Window um die relevante SD Koordinate,", 
       "die für die Berrechnung des SRP support genutzt wurde"),  
       
       paste0("nt up- and downstream Sequenze Window um die relevante SD",  
              "Koordinate mit der eingebauten Variation"),                
       
       paste0("nt up- and downstream Sequenze Window um die relevante SA", 
       "Koordinate "),                                            
       
       paste0("nt up- and downstream Sequenze Window um die relevante SA", 
              "Koordinate mit der eingebauten Variation"),              
       
       "Warum könnte diese Variation interessant sein?",     
       
       "neue Variation_ID" )
  
  docuDF <- data.frame(column=columns, documentation=documentation)
  
  return(docuDF)
}
