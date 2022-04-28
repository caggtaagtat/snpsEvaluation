createOutputDocu <- function(){
  
  columns <- c("gene",
               "c.DNA",
               "category",
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
               "SRE_support_alt_overlapping_SD",
               "MaxEnt alt overlapping SAs",
               "MaxEnt diff overlapping SAs",
               "MaxEnt diff overlapping SA vs annotated",
               "SRE_support_alt_overlapping_SA",
               "Maxent alt overlapping SD GC sites",
               "Maxent diff overlapping SD GC sites",
               "SRE_support_alt_overlapping_GCsite",
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
               "category details",
               "Next_SA_coordinate"=0,
               "Next_SD_coordinate"=0,
               "ID")
  
  documentation <- c(
       "gene name",  
       
       "Annotation der Variation",
       
       "Potentielle Art der Spleiss-Hinderung",
       
       "Differenz im HZEI-integral durch die variation",
       
       "Aenderungen im MaxEnt score des annotierten Akzeptor",
       
       "Ist das AG Dinucleotid des annotierten Akzeptor betroffen?",
                                                                                                                                                                                                                                                                                                               
       paste0("Aenderungen im HBond score des annotierten Donors. Ist 9999,", 
              "falls das GT des Donors durch die Variation gestoert wird"),                                                                                                                                                       
       
       "Ist das GT dinucleotid des annotierten Donors betroffen?",  
       
       "Entfernung der variation zum naechsten annotierten SD", 
       
       "Entfernung der Variation zum naechsten annotierten SA", 
       
       "Ist die Variation exonisch oder intronisch im Transkript",
       
       paste0("HZEI-Integral downstream minus HZEI-Integral upstream des", 
       "naechstgelegenen SA in Reference-Sequzenz"),   
       
       paste0("HZEI-Integral downstream minus HZEI-Integral upstream des", 
              "naechstgelegenen SA in Reference-Sequzenz mit der Variation"),                                                                                                                                                                                                                     
       
       paste0("Differenz im HZEI-Integral Gewicht des SA mit und ohne",  
       "Variation"),    
       
       paste0("(Differenz im HZEI-Integral Gewicht des SA) geteilt durch", 
              "SA SRE support in %"),              
       
       paste0("HZEI-Integral upstream minus HZEI-Integral downstream des",  
       "naechstgelegenen SD in Reference-Sequzenz"),          
       
       paste0("HZEI-Integral upstream minus HZEI-Integral downstream des", 
              "naechstgelegenen SD mit der Variation"),     
       
       "Differenz im HZEI-Integral Gewicht des SD", 
       
       paste0("(Differenz im HZEI-Integral Gewicht des SD) geteilt durch", 
       "SD SRE support in %"),                         
       
       paste0("Maximaler HBS von ueberlappenden SD-Sequenzen in",  
              "Referenz-Sequenz mit eingebauter Variation"),    
       
       paste0("Differenz des maximalen HBS von ueberlappenden",  
       "SD-Sequenzen. Unterschied nach Einbauen der Variation",
              "in Referenz-Sequenz"),     
       
       paste0("Differenz des maximalen HBS von ueberlappenden",  
       "SD-Sequenzen zu dem HBS des naechstegelegenen Donors.",  
       "Unterschied nach Einbauen der Variation
       in Referenz-Sequenz"), 
       
       paste0("SRE support der ueberlappenden SD-site mit dem hoechsten",  
              "HBond score"), 
       
       paste0("Maximaler MaxEnt score von ueberlappenden SA-Sequenzen", 
       "in Referenz-Sequenz mit eingebauter Variation"),     
       
       paste0("Differenz des maximalen MaxEnt scores  von ueberlappenden", 
              "SA-Sequenzen. Unterschied nach Einbauen der Variation", 
               "in Referenz-Sequenz"), 
       
       paste0("Differenz des maximalen MaxEnt scores von ueberlappenden",  
              "SA-Sequenzen zu dem MaxEnt scores des naechstegelegenen",  
       "Akzeptors Unterschied nach Einbauen der Variation", 
              "in Referenz-Sequenz"),
       
       paste0("SRE support der ueberlappenden SA-site mit dem hoechsten",  
              "MaxEnt score"), 
       
       paste0("Maximaler MaxEnt score von ueberlappenden SD-Sequenzen", 
       "mit einem GC-Dinukleotid auf SD Position +/+ in ", 
              "Referenz-Sequenz mit eingebauter Variation"),    
       
       paste0("Differenz des maximalen MaxEnt score  von ueberlappenden", 
       "SD-Sequenzen mit einem GC-Dinukleotid auf SD Position +/+.", 
              "Unterschied nach Einbauen der Variation in Referenz-Sequenz"),  
       
       paste0("SRE support der ueberlappenden GC-site mit dem hoechsten",  
              "HBond score"), 
       
       "HBS des annotierten SD in Referenzgenom",       
       
       paste0("HBS des annotierten SD mit der entsprechenden Variation", 
       "(-999 falls resultierende SD mit der Variation kein GT mehr hat)"), 
       
       "MaxEnt score des Akzeptors im Referenzgenom",  
       
       "MaxEnt score des Akzeptors mit der enstprechenden Variation", 
       
       paste0("Liegt der Variation upstream oder downstream zum naechsten", 
       "annotierten SD"),             
       
       paste0("Liegt der Variation upstream oder downstream zum naechsten", 
              "annotierten SA"),                                         
       
       "Laenge des naechst-gelegenen Exons",                   
       
       "Typ der Variation (Variation, DEL, DUP oder INS)", 
       
       paste0("Splice site Stärke der uebernaechsten splice site upstream der", 
       " Variation"),                                           
       
       paste0("Splice site Stärke der uebernaechsten splice site downstream der", 
              "Variation"),                                       
       
       "gene_strand",   
       
       "chromosome", 
       
       "ENSEMBL_transcriptID",
       
       "genome_position",   
       
       "alt",            
       
       "Exon_sequence" , 
       
       paste0("Umgebendes Sequenz-Window um die relevante SD Koordinate,", 
       "die fuer die Berrechnung des SRP support genutzt wurde"),  
       
       paste0("nt up- and downstream Sequenze Window um die relevante SD",  
              "Koordinate mit der eingebauten Variation"),                
       
       paste0("nt up- and downstream Sequenze Window um die relevante SA", 
       "Koordinate "),                                            
       
       paste0("nt up- and downstream Sequenze Window um die relevante SA", 
              "Koordinate mit der eingebauten Variation"),              
       
       "Warum koennte diese Variation interessant sein?",  
       
       "Koordinat des naechstegelegenen SA",
       
       "Koordinat des naechstegelegenen SD",
       
       "neue Variation_ID" )
  
  docuDF <- data.frame(column=columns, documentation=documentation)
  
  return(docuDF)
}
