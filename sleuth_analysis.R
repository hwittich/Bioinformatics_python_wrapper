library(sleuth)

args <- commandArgs(trailingOnly=TRUE)
if(args[1]=="True"){
  output_dir<-"test_outputs"
} else {
  output_dir<-"miniProject_Henry_Wittich"
}

#read in the table describing samples and kallisto output
stab <- read.table(output_dir+"/kallisto/sample_table.tsv",header=TRUE,stringsAsFactors = FALSE,sep="\t")
#sleuth object
so <- sleuth_prep(stab)

#fit a model comparing the two conditions
so <- sleuth_fit(so, ~days.post.infection, 'full')

#fit the reduced model to compare in the likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions
so <- sleuth_lrt(so,'reduced','full')

#Extracting results
library(dplyr)
sleuth_table <- sleuth_results(so,'reduced:full','lrt',show_all = FALSE)

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

#write the info for the significant SNPs to the output file
output <- sleuth_significant[,c(1,4,2,3)]
write.table(output, file=output_dir+"/miniProject.log",quote=FALSE,row.names=TRUE,sep="\t",append=TRUE)
