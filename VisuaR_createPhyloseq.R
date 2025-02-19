# # #http://web.stanford.edu/class/bios221/stamps/Lab_import.html
basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
package.list <- setdiff(package.list,basic.packages)
if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
rm(list=ls()); closeAllConnections(); gc()

VisuaRName <- ""
VisuaRPath <- ""

seq.tax <- data.table::fread(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/",VisuaRName,"_ASVbySample_abund_final.csv",sep=""))

n.samples <- length(readRDS(paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_Mgroups.rds",sep="")))

# # extract samples only
seq.table <- as.matrix(seq.tax[,3:(2+n.samples)]); rownames(seq.table) <- paste0("ASV",seq.tax$ASV)
tax.table <- as.matrix(seq.tax[,(3+n.samples):ncol(seq.tax)]); rownames(tax.table) <- paste0("ASV",seq.tax$ASV)
rm(seq.tax)

ASV <- phyloseq::otu_table(seq.table,taxa_are_rows = T)
TAX <- phyloseq::tax_table(tax.table)

context.data <- readRDS(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/",VisuaRName,"_DiversityIndices-t.rds"))

context.data <- phyloseq::sample_data(context.data)

physeq <- phyloseq::phyloseq(ASV,TAX,context.data)

saveRDS(physeq, paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_phyloseq.rds"))
