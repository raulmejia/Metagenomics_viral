######################
###
# Count table is a non headed, tab separated two column dataframe that contains:
### Reads or contigs in the first column
### Lane and sample in the second column
#### Extract of the count table ####
# Z48163.2 201001339_S1_L001
# Z48163.2 201001339_S1_L001
#
###  referencetable is a headed two column file that contains 
# In the first column the Sample&Lane
# In the second column the number of non-human reads
###  Extract of the reference library #####
# ID;numberofreads_after_bowtie2
# S1_L001;669188
# S1_L002;653369

######################
##
##########################
### Installing and or Loading the required packages 
##########################
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library("tidyverse")
}
if (!require("stringr")) {
  install.packages("stringr", ask =FALSE)
  library("stringr")
}
if (!require("ggplot2")) {
  install.packages("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("dplyr")) {
  install.packages("dplyr", ask =FALSE)
  library("dplyr")
}
#library(edgeR)
#library(data.table)

#######################
### Functions defined by the user 
#######################

############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()
# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]" )
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-t", "--counttabledf", type="character", 
                    help="input count table")
parser$add_argument( "-r", "--referencetable", type="character" , 
                     help="path to your reference table" )
parser$add_argument("-l", "--label", type="character", 
                    help="label to save your results")
parser$add_argument("-o", "--outputfolder", type="character", 
                    help="output file where you want to store your results")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## The program starts
#############################
mycountdf_path <-args$counttabledf
# mycountdf_path  <-"/media/rmejia/mountme88/Projects/Maja/Maja_PhD_virusproject/allfiles.csv"

referencetable_path <-args$referencetable
# referencetable_path <-"/media/rmejia/mountme88/Projects/Maja/Maja_PhD_virusproject/numberofreads_after_bowtie2.csv"

label <- args$label
# label <- "Contigs_from_filtered_reads"

save_folder <- args$outputfolder
#  save_folder <- "/media/rmejia/mountme88/Projects/Maja/Maja_PhD_virusproject/Results/Contigs_Normalized_by_Filtered_reads/"

####################
### The program starts
####################
mycountdf <- read.table(mycountdf_path, fill = T, stringsAsFactors = T, row.names = NULL) #skipNul = T
df_reference_library <- read.table( referencetable_path , header= T, sep = ";")

dir.create(save_folder, recursive = TRUE)
save_folder <- normalizePath(save_folder)

#allfilesmatrix <- as.matrix(allfiles)
colnames(mycountdf) <- c("sseqid", "Libsize_patient_lane")

# Splitting by patient
mycountdf$patient <- mycountdf$Libsize_patient_lane # Adding a column with the patient code.
mycountdf$patient <- gsub( "[[:graph:]]*_S","S" , mycountdf$patient)
mycountdf$patient <- gsub( "_[[:graph:]]*","" , mycountdf$patient)
mycountdf$patient <- as.factor(mycountdf$patient)

mycountdf_splitted_by_patient <- split( mycountdf ,mycountdf$patient )

list_counts_per_patient <- lapply( mycountdf_splitted_by_patient , function(x){ count(x, sseqid)} )


sseqid_countmatrix <- do.call( rbind, list_counts_per_patient ) # matrix of sseqidcounts
sseqid_countmatrix$samplenumber <- rownames(sseqid_countmatrix)
sseqid_countmatrix$samplenumber<- sub( "\\.[[:digit:]]*","", sseqid_countmatrix$samplenumber )
sseqid_countmatrix$libsizepersample <- rep("NA" , dim(sseqid_countmatrix)[1] ) 

# getting the reference library per patient 
df_reference_library$Patient <- sub( "_[[:graph:]]*","", df_reference_library$ID )
df_reference_library$Patient <- as.factor( df_reference_library$Patient )
splitted_per_patient_df_reference_library <- split( df_reference_library, df_reference_library$Patient)

list_reflibrary_per_patient <- lapply( splitted_per_patient_df_reference_library , function(x){sum(x[,2])} )
df_reference_per_patient <- data.frame(matrix(unlist(list_reflibrary_per_patient), nrow=length(list_reflibrary_per_patient), byrow=TRUE))
colnames( df_reference_per_patient) <- "lib_size"
df_reference_per_patient$patient <- names( list_reflibrary_per_patient )
df_reference_per_patient$numberofpatient <- df_reference_per_patient$patient
df_reference_per_patient$numberofpatient <- sub( "S", "", df_reference_per_patient$numberofpatient)
df_reference_per_patient$numberofpatient <- as.numeric( df_reference_per_patient$numberofpatient)


## Integrating both data sets
for( k in  1:dim(df_reference_per_patient)[1] ){
  pos <- which(sseqid_countmatrix[,"samplenumber"] %in%  df_reference_per_patient[k,"patient"]) 
  sseqid_countmatrix[pos,"libsizepersample"]<- df_reference_per_patient[k,"lib_size"]
}

# correction factor = x/mean()
sseqid_countmatrix$libsizepersample <- as.numeric(sseqid_countmatrix$libsizepersample)
sseqid_countmatrix$correctionfactor <- sseqid_countmatrix$libsizepersample/mean(sseqid_countmatrix$libsizepersample )

# Corrected
sseqid_countmatrix$n_over_itslibsize <- sseqid_countmatrix$n / sseqid_countmatrix$libsizepersample
sseqid_countmatrix$n_corrected_by_CF <- sseqid_countmatrix$n / sseqid_countmatrix$correctionfactor

# Saving matrix
matrix_CF_path <-paste0(save_folder,"/",label,"_Matrix_with_correction_factors_by_library_size.txt")
write.table(sseqid_countmatrix, file =matrix_CF_path , sep=";", col.names = TRUE, row.names = TRUE )

reference_per_patient_path <-paste0(save_folder,"/",label,"_Matrix_reference_per_patient.txt")
write.table( df_reference_per_patient, file =reference_per_patient_path , sep=";", col.names = TRUE, row.names = TRUE )


# Plotting
sseqid_countmatrix_splitted <- split( sseqid_countmatrix , sseqid_countmatrix$samplenumber )

# Ordering the plot object
sortme <-function(x){
  the_order <- order(x[,"n"] , decreasing=TRUE)
  x_sorted <- x[the_order,]
  return(x_sorted)
  }

sseqid_countmatrix_splitted_sorted <- lapply(sseqid_countmatrix_splitted , sortme)


#matrix <- sseqid_countmatrix_splitted[[1]]


  for ( k in 1:length( sseqid_countmatrix_splitted_sorted ) ) {
  pdf( file = paste0( save_folder,"/", label,"_",k,".pdf") , width = 7 , height = 7  )
    gg <-  ggplot( data= sseqid_countmatrix_splitted_sorted[[k]] , aes( x= reorder( sseqid, -n_corrected_by_CF),n_corrected_by_CF ))
    ggjingles <- gg + geom_col(aes(fill="n_corrected"))+
      theme(axis.text.x = element_text(angle = 65, hjust = 1), text = element_text(size=9)) +
      labs(x="Virus", y= "Counts corrected by libsize") +
      ggtitle(k)
  print(ggjingles)
  dev.off()
  }

for ( k in 1:length( sseqid_countmatrix_splitted_sorted ) ) {
  pdf( file = paste0( save_folder,"/", label,"_n_corrected_by_CF_",k,".pdf") , width = 7 , height = 7  )
    gg <-  ggplot( data= sseqid_countmatrix_splitted_sorted[[k]] , aes( x= reorder( sseqid, -n_corrected_by_CF),n_corrected_by_CF ))
    ggjingles <- gg + geom_col(aes(fill="n_corrected"))+
      theme(axis.text.x = element_text(angle = 65, hjust = 1), text = element_text(size=9)) +
      labs(x="Virus", y= "Counts corrected by libsize") +
      ggtitle(k)
  print(ggjingles)
  dev.off()
}
# 
# listofdf_2_matrix <- function( sseqid_countmatrix, column_to_extract_counts ){ 
#     
#   
#   }
# 
 Features_big_matrix <- unique(  as.character(  sseqid_countmatrix[,"sseqid"]  )  )
 Samples_big_matrix <- unique( as.character ( sseqid_countmatrix[,"samplenumber"] ) )

 BigMat <- matrix( # Building the matrix
   rep( NA, length( Features_big_matrix) * length( Samples_big_matrix ) ) ,
   nrow = length( Features_big_matrix)
   )
 rownames( BigMat) <- Features_big_matrix ; colnames( BigMat) <- Samples_big_matrix
 BigMat_corrected_by_libsize<- BigMat
 
 
 sseqid_countmatrix_splitted_sorted_sseqid_char <- lapply( sseqid_countmatrix_splitted_sorted , function(x){ x[,"sseqid"] <- as.character( x[,"sseqid"] ) ; return(x) } )
 sseqid_countmatrix_splitted_sorted_sseqid_char <- lapply( sseqid_countmatrix_splitted_sorted , function(x){ rownames(x) <- x[,"sseqid"]  ; return(x) } )
 
for( k in  names( sseqid_countmatrix_splitted_sorted_sseqid_char ) ) { 
  for ( w in sseqid_countmatrix_splitted_sorted_sseqid_char[[k]][,"sseqid"] ) {
    BigMat[ w, k] <- sseqid_countmatrix_splitted_sorted_sseqid_char[[k]][ w,"n"]
  }
}
 
 BigMat_path <-paste0(save_folder,"/",label,"_BigMat_UnNormalized_counts.tsv")
 write.table( BigMat , file = BigMat_path , sep="\t", col.names = TRUE, row.names = TRUE )
 
 
 for( k in  names( sseqid_countmatrix_splitted_sorted_sseqid_char ) ) { 
   for ( w in sseqid_countmatrix_splitted_sorted_sseqid_char[[k]][,"sseqid"] ) {
     BigMat_corrected_by_libsize[ w, k] <- sseqid_countmatrix_splitted_sorted_sseqid_char[[k]][ w,"n_corrected_by_CF"]
   }
 }
 BigMat_corrected_path <-paste0(save_folder,"/",label,"_BigMat_Corrected_by_library_size.tsv")
 write.table( BigMat_corrected_by_libsize , file = BigMat_corrected_path , sep="\t", col.names = TRUE, row.names = TRUE )
 

 BigMat_corrected_by_libsize[ is.na( BigMat_corrected_by_libsize) ] <- 0
 
pdf( file = paste0( save_folder,"/", label,"_BigMat_corrected_by_libsize_Abundant_patient_per_contig.pdf") , width = 7 , height = 7  )
heatmap(BigMat_corrected_by_libsize)   
dev.off()

pdf( file = paste0( save_folder,"/", label,"_BigMat_corrected_by_libsize_Abundant_ViralContigs_in_patients.pdf") , width = 7 , height = 7  )
heatmap( t (BigMat_corrected_by_libsize)  )  
dev.off()


# S8, S10, S11, S15 : all immunocompetent native kidneys but BK counts: 24, >60, 24, >25 respectively  high bitscores (high sequence similarities, so “wrong” hits unlikely)
# 12, S17, S18: native immunocompromised, BK count 26, >30, 24 could be undetected BK-Polyomanephropathy of native kidney
# S1, S2, S14: transplant but SV40- (SV40 is the Immunohistochemistry-stain for Polyomavirus): BK counts >20, >25, >20 respectively
# S3, S5, S7, S16 and S19: transplant, polyomavirus-positive controls
