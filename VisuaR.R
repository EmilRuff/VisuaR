# #=== VisuaR Info ========================================================================================================================================================================
# # Workflow by Emil Ruff & Isabella Hrabe de Angelis 05/2020
# # Copyright Emil Ruff
# # The authors acknowledge valuable input by Alban Ramette and Angelique Gobet
# # Please cite https://github.com/EmilRuff/VisuaR

# #=== 1. User INPUT =====================================================================================================================================================================

# #=== 1.1. Mandatory INPUT ===============================================================================================================================================================

VisuaRProjectName=''              # Give your analysis a (short) name

PathToVisuarOutput=file.path('')  # At this location the folder 'VisuaR' will be created containing subfolders with your analyses (named as provided in VisuaRProjectname), e.g. 'C:/Users/User1/Desktop'

# # Provide your dada2 output
PathToSeqtabNochim=file.path('')  # provide the location and name of your seqtab_nochim.rds, e.g. C:/Users/User1/Desktop/FancyProject_seqtab_nochim.rds
PathToTaxonomy=file.path('')      # provide the location and name of your taxonomy.rds, e.g. C:/Users/User1/Desktop/FancyProject_taxa_species.rds

# # Analysis specifications
KingdomOfInterest=''              # Which kingdom are you interested in? choose Bacteria/Archaea/Eukaryota
MinimumAllowedReadCount=2000      # MinimumAllowedReadcount depends on the overall quality of the run and targeted sequencing depth and is used to remove samples with low reads (failed samples, controls). All samples with less reads than this number will be excluded from the analysis. This is based on all reads derived from your dada2 analysis, no matter which kingdom they belong to. 

# #=== 1.2. Optional INPUT - Sample selection =============================================================================================================================================

CladeOfInterest=''                # Enter one or several clade(s) or keywords as found in your taxa.rds file, e.g. 'Desulfo|Meth'. Only keeps ASVs belonging to lineages containing the term 'Desulfo' and/or 'Meth', i.e. 'Methylococcus', 'Methylotenera'...
ExcludeClade=''                   # Enter one or several clade(s) or keywords as found in your taxa.rds file. e.g. 'Escherichia|Salmonella'. Excludes ASVs belonging to lineages containing the term 'Escherichia' and/or 'Salmonella', i.e. Escherichia coli.

KeepSamplesbyName=''              # Enter names of samples that you want to keep. Leave blank to keep all samples. For an exact match use e.g. '^SampleName1$|^SampleName2$'. The use of general terms is possible too, e.g. 'control|blank|enrichment|cont|cntrl|enr|extrctrl|Control|Blank|Enrichment|Cont|Cntrl|Enr|Extrctrl'.
excludeSamplesbyName=''           # Enter names of samples that you want to exclude. Leave blank to exclude no samples.

MinimumAllowedReadCount.Analysis=2000 # This field can be set to any number including 0 (zero) depending on which samples (based on read counts) should be included in the analysis. This happens after all ASVs and samples have been excluded as set in KeepSamplesbyName, excludeSamplesbyName, CladeOfInterest, ExcludeClade.

SaveWholeworkspace='Y'            # Set this variable to 'N' to save computational resources. R might run into memory issues if you have large datasets of e.g. 1000 samples. Set this variable to 'Y' to save all created tables in your workspace.

# #=== 1.3. Optional INPUT - Contextual data ==============================================================================================================================================

Metadata='N' # Do you have contextual data? If not write 'N' , if yes write 'Y' and select a grouping category (Grouping1) as found in your metadata header.

# # Useful notes on how to prepare your metadata file:
# # replace missing values by NA
# # replace spaces or special symbols by underscores
# # use unique row names, rownames have to be exactly the same as your sample names from your sequencing files
# # dates: use 01/01/2016 instead of 01/01/16
# # do not start row or column names with a number
# # use a .txt file

Metadatafilepath=file.path('')  # provide the file path to your contextual data. e.g. 'C:/Users/User1/FancyProject_metadata.txt' 
Grouping1=''
Grouping2=''                    # The samples within Grouping1 will be sorted using this continuous grouping parameter (e.g. concentration, depth profile) for visualization
VariableToExclude=''            # Enter a contextual data variable to exclude. E.g. '4C'
ColumnOfVariableToExclude=''    # Enter column name where the VariableToExclude can be found. E.g. 'Temperature'

# #=== 1.4. Specify VisuaR output =========================================================================================================================================================

NumberOfTOPClades=20      # Choose the number of most abundant clades that should be included in the relative abundance plot (bubbleplot). Remaining clades will be summed up and shown as 'Others'. A number between 10 and 30 is recommended.
NumberOfTOPClades.box=10  # Choose the number of most abundant clades that should be included in the relative abundance plot (boxplots). A number of 5 to 20 is recommended.

NI=10                     # number of iterations used for the calculation of alpha diversity indices (iterations will take long for big datasets, 10-50 are recommended based on number of samples analysed).

M.col=c('')               # You can provide a color vector, e.g. c('blue','orange') or c('#CC79A7','#999999'). The vector needs as many colors as unique variables found in 'Grouping1'. The colors will be assigned to the Grouping1 variables in alphabetical order. 

# #=== 2. Prologue =========================================================================================================================

# #=== 2.1. Creates file directories =========================================================================================================================
PathToVisuaRAnalysis=file.path(PathToVisuarOutput,'VisuaR',paste(VisuaRProjectName,sep='_')) # creates directories for the current VisuaR analysis
dir.create(file.path(PathToVisuarOutput,'VisuaR')) 
dir.create(PathToVisuaRAnalysis)
dir.create(file.path(PathToVisuaRAnalysis,'Alpha_Diversity'))
dir.create(file.path(PathToVisuaRAnalysis,'Beta_Diversity'))

# #=== 2.2. Creates log file =========================================================================================================================
cat('VISUAR Analysis',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\n\n1. User INPUT',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
sink(file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),append=TRUE,type='output')
cat('\n',sep='')
print(paste('AnalysisDate: ',Sys.Date()))
cat('VisuaRProjectName: ',VisuaRProjectName,'\n',sep='')
cat('Seqtab nochim originates from: ',PathToSeqtabNochim,'\n',sep='')
cat('Taxonomy originates from: ',PathToTaxonomy,'\n',sep='')
cat('Metadata originates from: ',Metadatafilepath,'\n\n',sep='')
cat('Samples belonging to the KingdomOfInterest: ',KingdomOfInterest,' will be kept for the further analysis.','\n',sep='')
if (MinimumAllowedReadCount!=''|MinimumAllowedReadCount!=0) {
  cat('Samples with less then ',MinimumAllowedReadCount,' reads (in all kingdoms, without filtering) will be excluded from the analysis.','\n',sep='')
} else {
  cat('No samples will be excluded due to low read counts (in all kingdoms, without filtering) will be excluded as no MinimumAllowedReadCount was set.','\n',sep='')
}
if (CladeOfInterest!='') {
  cat('Only ASVs belonging to the CladeOfInterest ',CladeOfInterest,' will be kept for further analysis','\n',sep='')
} else {
  cat('No CladeOfInterest was chosen.','\n',sep='')
}
if (ExcludeClade!='') {
  cat('ASVs belonging to ExcludeClade ',ExcludeClade,' will be excluded from the analysis.','\n',sep='')
} else {
  cat('No ExcludeClade was chosen.','\n',sep='')
}
if (KeepSamplesbyName!='') {
  cat('Samples will be kept based on: ','\n',sep='')
  cat('KeepSamplesbyName: ',KeepSamplesbyName,'\n',sep='')
} else {
  cat('All samples will be kept based on the sample name as no KeepSamplesbyName was chosen.','\n',sep='')
}
if (excludeSamplesbyName!='') {
  cat('Samples will be excluded based on: ','\n',sep='')
  cat('excludeSamplesbyName: ',excludeSamplesbyName,'\n',sep='')
} else {
  cat('No samples will be excluded based on the sample name as no excludeSamplesbyName was chosen.','\n',sep='')
}
if (exists('MinimumAllowedReadCount.Analysis')==F) {
  cat('No MinimumAllowedReadCount.Analysis was set or it was set to 0. No sample will be excluded based on too low read counts after subsetting.')
} else {
  cat('Samples with less than ',MinimumAllowedReadCount.Analysis,' reads (in the subsetted data) will be excluded from the analysis.','\n',sep='')
}
if (Metadata=='N') {
  cat('\nNo Metadata was provided.','\n',sep='')
} else {
  cat('\nMetadata was provided.','\n',sep='')
  cat('The metadata originated from ',Metadatafilepath,'\n',sep='')
  cat('The samples will be grouped based on the category ',Grouping1,' (Grouping1).','\n',sep='')
  if (Grouping2!='') {
    cat('The samples will be sorted inside Grouping1 based on the category ',Grouping2,' (Grouping2).','\n',sep='')
  } else {
    cat('No Grouping2 was provided. The samples will not be sorted inside the main Grouping1.','\n',sep='')
  }
  if (VariableToExclude!='') {
    cat('Samples belonging to the VariableToExclude: ',VariableToExclude,' as found in the contextual data column ',ColumnOfVariableToExclude,' will be excluded from the analysis','\n',sep='')
  } else {
    cat('No samples belonging to a particular category in your metadata sheet will be excluded as no VariableToExclude was chosen.','\n',sep='')
  }
}
sink()

cat('\n\n2. Prologue',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# #=== 2.3. Loads required packages =======================================================================================================================================================

cat('\n2.3. Load required packages',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nInformation on your used package Versions and R version can be found in the txt file VersionInformation.txt',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
sink(file=(file.path(PathToVisuaRAnalysis,'Session_Info.txt')),append=TRUE)
sink(stdout(),type='message')
if (!require('stringr')) {install.packages('stringr')}; library(stringr)
if (!require('plotrix')) {install.packages('plotrix')}; library(plotrix)
if (!require('plyr')) {install.packages('plyr')}; library(plyr)
if (!require('vegan')) {install.packages('vegan')}; library(vegan)
if (!require('reshape2')) {install.packages('reshape2')}; library(reshape2)
if (!require('ggplot2')) {install.packages('ggplot2')}; library(ggplot2)
if (!require('ggsignif')) {install.packages('ggsignif')}; library(ggsignif)
if (!require('EnvStats')) {install.packages('EnvStats')}; library(EnvStats)
if (!require('ggpubr')) {install.packages('ggpubr')}; library(ggpubr)
if (!require('ape')) {install.packages('ape')}; library(ape)
if (!require('UpSetR')) {install.packages('UpSetR')}; library(UpSetR)
if (!require('venn')) {install.packages('venn')}; library(venn)
sessionInfo()
closeAllConnections() # closes all currently open connections.


# #=== 2.4. Reads and converts dada2 sequence table and taxonomy to match format used by VisuaR ===============================================================================================
cat('\n\n2.4. Read and convert dada2 sequence table and taxonomy to match format used by VisuaR.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

if (file.exists(file.path(paste(str_remove(PathToTaxonomy, '.rds'),'_noNAs.rds',sep='')))) {
  M.taxo.noNA=readRDS(file=file.path(paste(str_remove(PathToTaxonomy, '.rds'),'_noNAs.rds',sep='')))
} else {
  M.taxo.noNA=readRDS(file=file.path(PathToTaxonomy))
}

M.taxo.noNA=as.data.frame(M.taxo.noNA)  # ASV by Taxonomy data frame
if (ncol(M.taxo.noNA)==6){              # adds a column to the taxonomy called 'Species' and fills the column with NAs if the taxonomy was classified on genus level. This is needed for VisuaR to work properly.
  sink(file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),append=TRUE,type='output')
  cat('\nYour taxonomy data is not on species level. A 7th column named Species, filled with NAs will be added.')
  sink()
  M.taxo.noNA$Species=rep(NA,nrow(M.taxo.noNA))
} else {
  sink(file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),append=TRUE,type='output')
  cat('\nYour taxonomy data is classfied on species level.')
  sink()
}

ncol.taxo.noNA=ncol(M.taxo.noNA) # Number of taxonomic levels. Should be 7.
nrow.taxo.noNA=nrow(M.taxo.noNA) # Number of rows in taxonomy file = number of observed ASVs

# #=== 2.5. Replaces NAs in taxonomy with 'unc'and concatenates the taxonomy ===============================================================================================================

if (!file.exists(file.path(paste(str_remove(PathToTaxonomy, '.rds'),'_noNAs.rds',sep='')))) { # Replaces all 'NA' by 'unc' and concatenates taxonomic levels using '_' as separator
  UnclassifiedTerm='unc'                                                                    
  M.taxo.noNA=as.matrix(M.taxo.noNA)
  M.taxo.noNA[is.na(M.taxo.noNA)]=UnclassifiedTerm 
  for(k in 1:nrow(M.taxo.noNA)){
    M.taxo.noNA[k,1]=(as.character(M.taxo.noNA[k,1]))
    M.taxo.noNA[k,2]=paste(as.character(M.taxo.noNA[k,1]),as.character(M.taxo.noNA[k,2]),sep='_')
    M.taxo.noNA[k,3]=paste(as.character(M.taxo.noNA[k,2]),as.character(M.taxo.noNA[k,3]),sep='_')
    M.taxo.noNA[k,4]=paste(as.character(M.taxo.noNA[k,3]),as.character(M.taxo.noNA[k,4]),sep='_')
    M.taxo.noNA[k,5]=paste(as.character(M.taxo.noNA[k,4]),as.character(M.taxo.noNA[k,5]),sep='_')
    M.taxo.noNA[k,6]=paste(as.character(M.taxo.noNA[k,5]),as.character(M.taxo.noNA[k,6]),sep='_')
    M.taxo.noNA[k,7]=paste(as.character(M.taxo.noNA[k,6]),as.character(M.taxo.noNA[k,7]),sep='_')
  }
  M.taxo.noNA=as.data.frame(M.taxo.noNA,stringsAsFactors=FALSE)
  saveRDS(M.taxo.noNA,file.path(file.path(paste(str_remove(PathToTaxonomy, '.rds'),'_noNAs.rds',sep=''))))
}

# #=== 2.6. Prepares Seqtab_nochim =========================================================================================================================================================

M.seqtab.nochim=readRDS(file.path(PathToSeqtabNochim)) 
ncol.seqtab.nochim=ncol(M.seqtab.nochim) # Number of columns in Seqtab_nochim = number of observed ASVs

ASV.names=rownames(M.taxo.noNA)
if (identical(rownames(M.taxo.noNA),colnames(M.seqtab.nochim))){                # checks whether the ASVs in the taxonomy and in the seqtab_nochim are identical
  M.taxo.noNA=data.frame(lapply(M.taxo.noNA, function(x){gsub('[()]','-',x)}))  # replaces all brackets in the taxonomy file with a '-'. While this happens the ASVs obtain a running number instead of the ASV Sequence.
  colnames(M.seqtab.nochim)=rownames(M.taxo.noNA)                               # replaces ASV identifiers in M.seqtab.nochim by running number as well
  M.seqtab.nochim=as.matrix(as.data.frame(M.seqtab.nochim)) 
} else {
  cat("\nError: You're Taxonomy file does not match your Seqtab.\n",file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\t',append=TRUE)
  stop("You're Taxonomy file does not match your Seqtab.")
}

ASV.mapfile=cbind(ASV.names,M.taxo.noNA)                                            # prints a summary to archive which ASV identifier (running number) belongs to which ASV sequence
ASV.mapfile=cbind(row.names(ASV.mapfile),data.frame(ASV.mapfile,row.names = NULL))  # reformats the table to obtain excel readable output file
colnames(ASV.mapfile)[1]='ASV'
colnames(ASV.mapfile)[2]='ASV sequence'
write.table(ASV.mapfile,file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_ASV_summary.txt',sep='')),sep='\t',col.names = NA)

if (SaveWholeworkspace=='N') { # removes file for efficient use of resources
  rm(ASV.mapfile,ASV.names) 
}

# #=== 2.7. Merges seqtab_nochim and taxonomy ==============================================================================================================================================

M.seq.tax=cbind.data.frame(t(M.seqtab.nochim),M.taxo.noNA)  # creates dataframe with all samples, ASVs and taxonomy. 
if (SaveWholeworkspace=='N') {rm(M.taxo.noNA)}

ncol.seq.tax=ncol(M.seq.tax)                                # Number of samples + number of taxonomy classes (7)

M.reads=rowSums(M.seqtab.nochim[,])                         # Number of reads per sample

M.ASV.reads=colSums(M.seqtab.nochim[,]) # Number of reads per ASV
M.ASV.reads.df=as.data.frame(M.ASV.reads)
M.ASV.reads.df=cbind('ASV'=rownames(M.ASV.reads.df),M.ASV.reads.df)
rownames(M.ASV.reads.df)=NULL 
colnames(M.ASV.reads.df)[2]='Reads'
write.table(M.ASV.reads.df,file.path(PathToVisuaRAnalysis,'Alpha_Diversity',paste(VisuaRProjectName,'_ReadsperASV_Original.txt',sep='')),sep='\t',col.names = NA)
if (SaveWholeworkspace=='N') {rm(M.ASV.reads.df)}

sink(file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),append=TRUE,type='output')
cat('\n\nGeneral Information on sequenced samples.\n',sep=' ')
cat('You started with a total of ',nrow(M.seqtab.nochim),' samples.\n',sep='')
cat("You're original samples contain in total",sum(M.reads),'ASV counts (reads).\n',sep=' ')
cat("You're original samples contain on average",round(mean(M.reads),1),'ASV counts, with a minimum of',min(M.reads),'ASV counts and a maximum of',max(M.reads),'ASV counts.\n',sep=' ')
cat("You're original ASVs occur on average",round(mean(M.ASV.reads),1),'times, with a minimum of',min(M.ASV.reads),'and a maximum of',max(M.ASV.reads),'times.\n',sep=' ')
sink()
if (SaveWholeworkspace=='N') {rm(M.ASV.reads,M.reads)}

M.differentASVsperSample=rowSums(M.seqtab.nochim[,]!= 0) # creates a table with the number of different ASVs per sample
sink(file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),append=TRUE,type='output')
cat("You're original samples contain on average",round(mean(M.differentASVsperSample),1),'different ASVs (ASV richness), with a minimum of',min(M.differentASVsperSample),'and a maximum of',max(M.differentASVsperSample),'.\n',sep=' ')
sink()

M.reads.df=as.data.frame(M.reads) 
M.differentASVsperSample=as.data.frame(M.differentASVsperSample)
M.differentASVsperSample=cbind(row.names(M.seqtab.nochim),data.frame(M.reads.df,row.names = NULL),data.frame(M.differentASVsperSample,row.names = NULL))
colnames(M.differentASVsperSample)[1]='Sample Name'
colnames(M.differentASVsperSample)[2]='Total Reads'
colnames(M.differentASVsperSample)[3]='Observed ASVs - Richness'
write.table(M.differentASVsperSample,file.path(PathToVisuaRAnalysis,'Alpha_Diversity',paste(VisuaRProjectName,'_ReadsandASVsperSample_Original.txt',sep='')),sep='\t',col.names = NA )
if (SaveWholeworkspace=='N') {rm(M.differentASVsperSample)}

# #=== 2.8. Excludes samples based on the MinimumAllowedReadCount ==========================================================================================================================
cat('\n2.8. Excludes samples based on read counts lower than ',MinimumAllowedReadCount,' (MinimumAllowedReadCount).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (MinimumAllowedReadCount==''|MinimumAllowedReadCount==0) {
  cat('\nNo samples were excluded based on their read counts as no MinimumAllowedReadCount was chosen.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  M.seq.tax.subset=M.seq.tax
} else {
  M.seq.tax.MinimumAllowedReadCount.samples.excluded=as.data.frame(M.seq.tax[,-which(numcolwise(sum)(M.seq.tax)>=MinimumAllowedReadCount)]) 
  if (ncol(M.seq.tax.MinimumAllowedReadCount.samples.excluded)==ncol.taxo.noNA) { # this is TRUE if all samples have read counts higher than MinimumAllowedReadCount
    cat('\nNo samples were excluded based on low read counts as all samples had more than ',MinimumAllowedReadCount,' reads.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    M.seq.tax.subset=M.seq.tax
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.samples.excluded)}
  } else {
    M.seq.tax.subset=M.seq.tax[,-which(numcolwise(sum)(M.seq.tax)<MinimumAllowedReadCount)] # only keeps samples with read counts higher than or equal to MinimumAllowedReadCount
    cat('\nThe following ',(ncol(M.seq.tax.MinimumAllowedReadCount.samples.excluded)-ncol.taxo.noNA),' samples were exclulded based on low read counts:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat(colnames(M.seq.tax.MinimumAllowedReadCount.samples.excluded[1:((ncol(M.seq.tax.MinimumAllowedReadCount.samples.excluded))-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=TRUE)
    M.seq.tax.MinimumAllowedReadCount.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # makes a df of ASVs which were removed from the dataset.
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.samples.excluded)}
    if (nrow(M.seq.tax.MinimumAllowedReadCount.ASVs.excluded)!=0) {
      M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] 
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.ASVs.excluded),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=',',append=T)
      if(SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.ASVs.excluded)}
    } else {
      cat('No ASVs have been completely excluded.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      if(SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.ASVs.excluded)}
    }
  }
}
cat((ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.9. Excludes ASVS not belonging to the KingdomofInterest ===========================================================================================================================

cat('\n\n2.9. Excludes ASVs not belonging to the chosen kingdom ',KingdomOfInterest,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (KingdomOfInterest=='') {
  cat('\nError: You did not choose a kingdom of interest.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\t',append=TRUE)
  stop('You did not choose a kingdom of interest.')
} else {
  M.seq.tax.kingdom.ASV.excluded.org=M.seq.tax[grep(KingdomOfInterest,M.seq.tax[,(ncol.seq.tax-6)],invert=T),]                                      # creates an ASV by Sample&Taxonomy df with ASVs not belonging to the KingdomOfInterest (based on the original dataset (M.seq.tax))
  M.seq.tax.kingdom.ASV.excluded=as.data.frame(M.seq.tax.subset[grep(KingdomOfInterest,M.seq.tax.subset[,(ncol(M.seq.tax.subset)-6)],invert=T),])   # creates an ASV by Sample&Taxonomy df with ASVs not belonging to the KingdomOfInterest (based on the subsetted M.seq.tax.subset)
  M.seq.tax.kingdom.ASV.kept.org=M.seq.tax[grep(KingdomOfInterest,M.seq.tax[,(ncol.seq.tax-6)],invert=F),]                                          # creates an ASV by Sample&Taxonomy df with ASVs and Samples belonging to the KingdomOfInterst (based on the original dataset (M.seq.tax))
  M.seq.tax.subset=M.seq.tax.subset[grep(KingdomOfInterest,M.seq.tax.subset[,(ncol(M.seq.tax.subset)-6)],invert=F),]                                # Further Subsets M.seq.tax.subset, such that only ASVs belonging to the KingdomOfInterest are present (based on previous subsetting)
  cat('\nThe following ',nrow(M.seq.tax.kingdom.ASV.excluded.org),' out of ',nrow(M.seq.tax),' ASVs have been excluded because they do not belong to the chosen kingdom ',KingdomOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat(rownames(M.seq.tax.kingdom.ASV.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
  cat('\n',nrow(M.seq.tax.subset),' out of ',nrow(M.seq.tax),' ASVs remain belonging to the kingdom ',KingdomOfInterest,'.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  M.seq.tax.kingdom.samples.excluded.org=M.seq.tax.kingdom.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.kingdom.ASV.kept.org)!=0)]                # creates an ASV by Sample&Taxonomy df of samples which no longer have any ASV counts belonging to the chosen KingdomOfInterest (based on the original dataset (M.seq.tax))
  M.seq.tax.kingdom.samples.excluded=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)!=0)]                                                # creates an ASV by Sample&Taxonomy df of samples which no longer have any ASV counts belonging to the chosen KingdomOfInterest (based on the subsetted M.seq.tax.subset)
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.kingdom.ASV.excluded,M.seq.tax.kingdom.ASV.excluded.org)}
  if (ncol(M.seq.tax.kingdom.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.kingdom.samples.excluded.org)!=ncol.taxo.noNA) {                      # only happens if at least one sample was excluded in both dfs (M.seq.tax and M.seq.tax.subset)
    M.seq.tax.kingdom.ASV.kept.org=M.seq.tax.kingdom.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.kingdom.ASV.kept.org)==0)]                      # excludes samples which no longer have any ASV counts in original df
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)]                                                                # excludes samples which no longer have any ASV counts in subsetted df
    cat('\nThe following ',(ncol(M.seq.tax.kingdom.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen kingdom ',KingdomOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.kingdom.samples.excluded.org[c(1:(ncol(M.seq.tax.kingdom.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append = T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.kingdom.ASV.kept.org,M.seq.tax.kingdom.samples.excluded,M.seq.tax.kingdom.samples.excluded.org)}
  } else if (ncol(M.seq.tax.kingdom.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.kingdom.samples.excluded.org)==ncol.taxo.noNA) {               # only happens if at least one sample was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)]                                                                # excludes samples which no longer have any ASV counts in subsetted df
    cat('\nThe following ',(ncol(M.seq.tax.kingdom.samples.excluded)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen kingdom ',KingdomOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.kingdom.samples.excluded[c(1:(ncol(M.seq.tax.kingdom.samples.excluded)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append = T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.kingdom.ASV.kept.org,M.seq.tax.kingdom.samples.excluded,M.seq.tax.kingdom.samples.excluded.org)}
  } else if (ncol(M.seq.tax.kingdom.samples.excluded)==ncol.taxo.noNA&ncol(M.seq.tax.kingdom.samples.excluded.org)!=ncol.taxo.noNA) {               # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.kingdom.ASV.kept.org=M.seq.tax.kingdom.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.kingdom.ASV.kept.org)==0)]                      # excludes samples which no longer have any ASV counts in original df
    cat('\nThe following ',(ncol(M.seq.tax.kingdom.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen kingdom ',KingdomOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.kingdom.samples.excluded.org[c(1:(ncol(M.seq.tax.kingdom.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append = T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.kingdom.ASV.kept.org,M.seq.tax.kingdom.samples.excluded,M.seq.tax.kingdom.samples.excluded.org)}
  } else if (ncol(M.seq.tax.kingdom.samples.excluded)==ncol.taxo.noNA&ncol(M.seq.tax.kingdom.samples.excluded.org)==ncol.taxo.noNA) {               # only happens if no sample was excluded
    cat('No samples were excluded due to no observed ASVs in the chosen kingdom.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.kingdom.ASV.kept.org,M.seq.tax.kingdom.samples.excluded,M.seq.tax.kingdom.samples.excluded.org)}
  }
}
cat((ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.10. Keeps particular ASVs based on the clade they belong to ========================================================================================================================
cat('\n\n2.10. Keeps particular ASVs based on the clade they belong to',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (CladeOfInterest=='') {
  cat('\nNo Clade of Interest was chosen.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append = T)
} else {
  M.seq.tax.clade.ASV.excluded.org=M.seq.tax[grep(CladeOfInterest,paste(M.seq.tax[,ncol(M.seq.tax)],M.seq.tax[,ncol(M.seq.tax)-1],M.seq.tax[,ncol(M.seq.tax)-2],M.seq.tax[,ncol(M.seq.tax)-3],M.seq.tax[,ncol(M.seq.tax)-4],M.seq.tax[,ncol(M.seq.tax)-5],M.seq.tax[,ncol(M.seq.tax)-6]),invert=T),]   # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the CladeOfInterest (based on the original dataset (M.seq.tax)) # would also work because all classifiers are in the Species column:M.seq.tax.clade.excluded.org=M.seq.tax[grep(CladeOfInterest,M.seq.tax[,ncol(M.seq.tax)],invert=T),]
  M.seq.tax.clade.ASV.excluded=M.seq.tax.subset[grep(CladeOfInterest,paste(M.seq.tax.subset[,ncol(M.seq.tax.subset)],M.seq.tax.subset[,ncol(M.seq.tax.subset)-1],M.seq.tax.subset[,ncol(M.seq.tax.subset)-2],M.seq.tax.subset[,ncol(M.seq.tax.subset)-3],M.seq.tax.subset[,ncol(M.seq.tax.subset)-4],M.seq.tax.subset[,ncol(M.seq.tax.subset)-5],M.seq.tax.subset[,ncol(M.seq.tax.subset)-6]),invert=T),] # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the CladeOfInterest (based on the subsetted M.seq.tax.subset)
  M.seq.tax.clade.ASV.kept.org=M.seq.tax[grep(CladeOfInterest,paste(M.seq.tax[,ncol(M.seq.tax)],M.seq.tax[,ncol(M.seq.tax)-1],M.seq.tax[,ncol(M.seq.tax)-2],M.seq.tax[,ncol(M.seq.tax)-3],M.seq.tax[,ncol(M.seq.tax)-4],M.seq.tax[,ncol(M.seq.tax)-5],M.seq.tax[,ncol(M.seq.tax)-6]),invert=F),] # creates an ASV by Sample&Taxonomy df with only those samples belonging to the CladeOfInterest (based on the original dataset (M.seq.tax))
  M.seq.tax.subset=M.seq.tax.subset[grep(CladeOfInterest,paste(M.seq.tax.subset[,ncol(M.seq.tax.subset)],M.seq.tax.subset[,ncol(M.seq.tax.subset)-1],M.seq.tax.subset[,ncol(M.seq.tax.subset)-2],M.seq.tax.subset[,ncol(M.seq.tax.subset)-3],M.seq.tax.subset[,ncol(M.seq.tax.subset)-4],M.seq.tax.subset[,ncol(M.seq.tax.subset)-5],M.seq.tax.subset[,ncol(M.seq.tax.subset)-6]),invert=F),] # Further Subsets M.seq.tax.subset, such that only ASVs belonging to the CladeOfInterest are present (based on previous subsetting)
  cat('\n',nrow(M.seq.tax.clade.ASV.excluded.org),' out of ',nrow(M.seq.tax),' ASVs have been excluded because they do not belong to the chosen clade.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat('\n',nrow(M.seq.tax.subset),' out of ',nrow(M.seq.tax),' ASVs remain belonging to the clade ',CladeOfInterest,' or where this expression occurs in the taxonomy.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  M.seq.tax.clade.samples.excluded.org=M.seq.tax.clade.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade.ASV.kept.org)!=0)] # creates an ASV by Sample&Taxonomy df of all samples which no longer have any ASV counts belonging to the chosen CladeOfInterest (based on the original dataset (M.seq.tax))
  M.seq.tax.clade.samples.excluded=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)!=0)] # creates an ASV by Sample&Taxonomy df of all samples which no longer have any ASV counts belonging to the chosen CladeOfInterest (based on the subsetted M.seq.tax.subset)
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade.ASV.excluded.org,M.seq.tax.clade.ASV.excluded)}
  if (ncol(M.seq.tax.clade.samples.excluded)!=ncol.taxo.noNA&(ncol(M.seq.tax.clade.samples.excluded.org)!=ncol.taxo.noNA)) { # only happens if at least one sample was excluded in both dfs
    M.seq.tax.clade.ASV.kept.org=M.seq.tax.clade.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade.ASV.kept.org)==0)] # excludes samples which no longer have any ASV counts
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen CladeOfInterest ',CladeOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade.samples.excluded.org[c(1:(ncol(M.seq.tax.clade.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade.ASV.kept.org,M.seq.tax.clade.samples.excluded.org,M.seq.tax.clade.samples.excluded)}
  } else if (ncol(M.seq.tax.clade.samples.excluded)!=ncol.taxo.noNA&(ncol(M.seq.tax.clade.samples.excluded.org)==ncol.taxo.noNA)) {# only happens if at least one sample was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade.samples.excluded)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen CladeOfInterest ',CladeOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade.samples.excluded[c(1:(ncol(M.seq.tax.clade.samples.excluded)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade.ASV.kept.org,M.seq.tax.clade.samples.excluded.org,M.seq.tax.clade.samples.excluded)}
  } else if (ncol(M.seq.tax.clade.samples.excluded)==ncol.taxo.noNA&(ncol(M.seq.tax.clade.samples.excluded.org)!=ncol.taxo.noNA)) { # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.clade.ASV.kept.org=M.seq.tax.clade.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade.ASV.kept.org)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the chosen CladeOfInterest ',CladeOfInterest,'.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade.samples.excluded.org[c(1:(ncol(M.seq.tax.clade.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade.ASV.kept.org,M.seq.tax.clade.samples.excluded.org,M.seq.tax.clade.samples.excluded)}
  } else if (ncol(M.seq.tax.clade.samples.excluded)==ncol.taxo.noNA&(ncol(M.seq.tax.clade.samples.excluded.org)==ncol.taxo.noNA)) { # only happens if at no sample was excluded
    cat('\nNo samples were excluded due to no observed ASVs in the chosen clade.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade.ASV.kept.org,M.seq.tax.clade.samples.excluded.org,M.seq.tax.clade.samples.excluded)}
  }
}
cat('\n',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.11. Excludes particular ASVs based on the clade they belong to =====================================================================================================================

cat('\n\n2.11. Excludes particular ASVs based on the clade they belong to.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (ExcludeClade=='') {
  cat('\nNo Clade to exclude was chosen.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append = T)
} else {
  M.seq.tax.clade2.ASV.excluded.org=M.seq.tax[grep(ExcludeClade,paste(M.seq.tax[,ncol(M.seq.tax)],M.seq.tax[,ncol(M.seq.tax)-1],M.seq.tax[,ncol(M.seq.tax)-2],M.seq.tax[,ncol(M.seq.tax)-3],M.seq.tax[,ncol(M.seq.tax)-4],M.seq.tax[,ncol(M.seq.tax)-5],M.seq.tax[,ncol(M.seq.tax)-6]),invert=F),]  # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'ExcludeClade' (based on the original dataset (M.seq.tax))
  M.seq.tax.clade2.ASV.excluded=M.seq.tax.subset[grep(ExcludeClade,paste(M.seq.tax.subset[,ncol(M.seq.tax.subset)],M.seq.tax.subset[,ncol(M.seq.tax.subset)-1],M.seq.tax.subset[,ncol(M.seq.tax.subset)-2],M.seq.tax.subset[,ncol(M.seq.tax.subset)-3],M.seq.tax.subset[,ncol(M.seq.tax.subset)-4],M.seq.tax.subset[,ncol(M.seq.tax.subset)-5],M.seq.tax.subset[,ncol(M.seq.tax.subset)-6]),invert=F),] # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'ExcludeClade' (based on the subsetted M.seq.tax.subset)
  M.seq.tax.clade2.ASV.kept.org=M.seq.tax[grep(ExcludeClade,paste(M.seq.tax[,ncol(M.seq.tax)],M.seq.tax[,ncol(M.seq.tax)-1],M.seq.tax[,ncol(M.seq.tax)-2],M.seq.tax[,ncol(M.seq.tax)-3],M.seq.tax[,ncol(M.seq.tax)-4],M.seq.tax[,ncol(M.seq.tax)-5],M.seq.tax[,ncol(M.seq.tax)-6]),invert=T),]  # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the 'Excluded clade' (based on the original dataset (M.seq.tax))
  M.seq.tax.subset=M.seq.tax.subset[grep(ExcludeClade,paste(M.seq.tax.subset[,ncol(M.seq.tax.subset)],M.seq.tax.subset[,ncol(M.seq.tax.subset)-1],M.seq.tax.subset[,ncol(M.seq.tax.subset)-2],M.seq.tax.subset[,ncol(M.seq.tax.subset)-3],M.seq.tax.subset[,ncol(M.seq.tax.subset)-4],M.seq.tax.subset[,ncol(M.seq.tax.subset)-5],M.seq.tax.subset[,ncol(M.seq.tax.subset)-6]),invert=T),] # Further Subsets M.seq.tax.subset, such that only ASVs not belonging to the 'ExcludeClade' are present (based on previous subsetting)
  cat('\n',nrow(M.seq.tax.clade2.ASV.excluded.org),' out of ',nrow(M.seq.tax),' ASVs have been excluded because they do belong to the clade ',ExcludeClade,'.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat('\n',nrow(M.seq.tax.subset),' out of ',nrow(M.seq.tax),' ASVs remain not belonging to the clade ',ExcludeClade,'.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  M.seq.tax.clade2.samples.excluded.org=M.seq.tax.clade2.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade2.ASV.kept.org)!=0)] # creates an ASV by Sample&Taxonomy df of all samples which no longer have any ASV counts (based on the original dataset (M.seq.tax))
  M.seq.tax.clade2.samples.excluded=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)!=0)] # creates an ASV by Sample&Taxonomy df of all samples which no longer have any ASV counts (based on the subsetted M.seq.tax.subset)
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade2.ASV.excluded.org,M.seq.tax.clade2.ASV.excluded)}
  if (ncol(M.seq.tax.clade2.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.clade2.samples.excluded.org)!=ncol.taxo.noNA) { # only happens if at least one sample was excluded in both dfs
    M.seq.tax.clade2.ASV.kept.org=M.seq.tax.clade2.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade2.ASV.kept.org)==0)] # excludes samples which no longer have any ASV counts
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade2.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the remaining clades.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade2.samples.excluded.org[c(1:(ncol(M.seq.tax.clade2.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade2.ASV.kept.org,M.seq.tax.clade2.samples.excluded.org,M.seq.tax.clade2.samples.excluded)}
  } else if (ncol(M.seq.tax.clade2.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.clade2.samples.excluded.org)==ncol.taxo.noNA) {  # only happens if at least one sample was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade2.samples.excluded)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the remaining clades.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade2.samples.excluded[c(1:(ncol(M.seq.tax.clade2.samples.excluded)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade2.ASV.kept.org,M.seq.tax.clade2.samples.excluded.org,M.seq.tax.clade2.samples.excluded)}
  } else if (ncol(M.seq.tax.clade2.samples.excluded)==ncol.taxo.noNA&ncol(M.seq.tax.clade2.samples.excluded.org)!=ncol.taxo.noNA) { # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.clade2.ASV.kept.org=M.seq.tax.clade2.ASV.kept.org[,-which(numcolwise(sum)(M.seq.tax.clade2.ASV.kept.org)==0)] # excludes samples which no longer have any ASV counts
    cat('\nThe following ',(ncol(M.seq.tax.clade2.samples.excluded.org)-ncol.taxo.noNA),' samples have been excluded due to no observed ASVs in the remaining clades.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.clade2.samples.excluded.org[c(1:(ncol(M.seq.tax.clade2.samples.excluded.org)-ncol.taxo.noNA))]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade2.ASV.kept.org,M.seq.tax.clade2.samples.excluded.org,M.seq.tax.clade2.samples.excluded)}
  } else if (ncol(M.seq.tax.clade2.samples.excluded)==ncol.taxo.noNA&ncol(M.seq.tax.clade2.samples.excluded.org)==ncol.taxo.noNA) { # only happens if no sample was excluded
    cat('\nNo samples were excluded due to no observed ASVs in the remaining clades.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.clade2.ASV.kept.org,M.seq.tax.clade2.samples.excluded.org,M.seq.tax.clade2.samples.excluded)}
  }
}
cat('\n',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.12. Keeps particular samples based on the sample name ==============================================================================================================================

cat('\n\n2.12. Keeps particular samples based on the sample name.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (KeepSamplesbyName=='') {
  cat('\nNo particular samples to keep were chosen.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append = T)
} else {
  M.seq.tax.keepsample.excluded.org=M.seq.tax[,grep(KeepSamplesbyName,colnames(M.seq.tax),invert=T)] # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the 'KeepSamplesbyName' (includes taxonomy, based on the original dataset (M.seq.tax))
  M.seq.tax.keepsample.excluded=as.data.frame(M.seq.tax.subset[,grep(KeepSamplesbyName,colnames(M.seq.tax.subset),invert=T)])  # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the 'KeepSamplesbyName' (includes taxonomy, based on the subsetted dataset (M.seq.tax.subset))
  M.seq.tax.keepsample.kept.org=cbind(M.seq.tax[,grep(KeepSamplesbyName,colnames(M.seq.tax),invert=F)],M.seq.tax[,(ncol(M.seq.tax)-ncol.taxo.noNA+1):ncol(M.seq.tax)]) # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'KeepSamplesByName' (based on the original dataset (M.seq.tax)). The taxonomy is deleted by the grep function and therefore added with cbind
  M.seq.tax.subset=cbind(M.seq.tax.subset[,grep(KeepSamplesbyName,colnames(M.seq.tax.subset),invert=F)],M.seq.tax.subset[,(ncol(M.seq.tax.subset)-ncol.taxo.noNA+1):ncol(M.seq.tax.subset)]) # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'KeepSamplesByName' (based on the subsetted dataset (M.seq.tax.subset)). The taxonomy is deleted by the grep function and therefore added with cbind
  cat('\n','Based on the expressions chosen to keep:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat(KeepSamplesbyName,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
  cat('The following ',ncol(M.seq.tax.keepsample.excluded.org)-ncol.taxo.noNA,' samples were excluded.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat(colnames(M.seq.tax.keepsample.excluded.org[,1:(ncol(M.seq.tax.keepsample.excluded.org)-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
  M.seq.tax.keepsample.ASVs.excluded.org=M.seq.tax.keepsample.kept.org[rowSums(M.seq.tax.keepsample.kept.org[,1:(ncol(M.seq.tax.keepsample.kept.org)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the original dataset (M.seq.tax))
  M.seq.tax.keepsample.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the subsetted dataset (M.seq.tax.subset))
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.keepsample.excluded.org,M.seq.tax.keepsample.excluded)}
  if (nrow(M.seq.tax.keepsample.ASVs.excluded)!=0&nrow(M.seq.tax.keepsample.ASVs.excluded.org)!=0) { # only happens if at least one ASV was excluded in both dfs
    M.seq.tax.keepsample.kept.org=M.seq.tax.keepsample.kept.org[rowSums(M.seq.tax.keepsample.kept.org[,1:(ncol(M.seq.tax.keepsample.kept.org)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('The following ',nrow(M.seq.tax.keepsample.ASVs.excluded.org),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.keepsample.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.keepsample.kept.org,M.seq.tax.keepsample.ASVs.excluded.org,M.seq.tax.keepsample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.keepsample.ASVs.excluded)!=0&nrow(M.seq.tax.keepsample.ASVs.excluded.org)==0) {  # only happens if at least one ASv was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('The following ',nrow(M.seq.tax.keepsample.ASVs.excluded),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.keepsample.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.keepsample.kept.org,M.seq.tax.keepsample.ASVs.excluded.org,M.seq.tax.keepsample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.keepsample.ASVs.excluded)==0&nrow(M.seq.tax.keepsample.ASVs.excluded.org)!=0) { # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.keepsample.kept.org=M.seq.tax.keepsample.kept.org[rowSums(M.seq.tax.keepsample.kept.org[,1:(ncol(M.seq.tax.keepsample.kept.org)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('The following ',nrow(M.seq.tax.keepsample.ASVs.excluded.org),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.keepsample.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.keepsample.kept.org,M.seq.tax.keepsample.ASVs.excluded.org,M.seq.tax.keepsample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.keepsample.ASVs.excluded)==0&nrow(M.seq.tax.keepsample.ASVs.excluded.org)==0) { # only happens if no ASV was excluded
    cat('\nNo ASVs were excluded due to no observed ASVs in the remaining samples.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.keepsample.kept.org,M.seq.tax.keepsample.ASVs.excluded.org,M.seq.tax.keepsample.ASVs.excluded)}
  }
}
cat('\n',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.13. Excludes particular samples based on the sample name ===========================================================================================================================

cat('\n\n2.13. Excludes particular samples based on the sample name.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

if (excludeSamplesbyName=='') {
  cat('\nNo particular samples to exclude were chosen.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append = T)
} else {
  M.seq.tax.excludesample.excluded.org=M.seq.tax[,grep(excludeSamplesbyName,colnames(M.seq.tax),invert=F)]                                  # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'ExcludeSamplesbyName' (based on the original dataset (M.seq.tax))
  M.seq.tax.excludesample.excluded=as.data.frame(M.seq.tax.subset[,grep(excludeSamplesbyName,colnames(M.seq.tax.subset),invert=F)])         # creates an ASV by Sample&Taxonomy df with only those samples belonging to the 'ExcludeSamplebyName' (based on the subsetted dataset (M.seq.tax.subset))
  M.seq.tax.excludesample.kept.org=M.seq.tax[,grep(excludeSamplesbyName,colnames(M.seq.tax),invert=T)]                                      # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the 'ExcludeSamplebyName' (based on the original dataset (M.seq.tax)). 
  M.seq.tax.subset=M.seq.tax.subset[,grep(excludeSamplesbyName,colnames(M.seq.tax.subset),invert=T)]                                        # creates an ASV by Sample&Taxonomy df with only those samples not belonging to the 'ExcludeSamplebyName' (based on the subsetted dataset (M.seq.tax.subset)).
  cat('\n','Based on the expressions chosen to exclude:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat(excludeSamplesbyName,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
  cat('The following ',ncol(M.seq.tax.excludesample.excluded.org),' samples were excluded.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
  cat(colnames(M.seq.tax.excludesample.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
  M.seq.tax.excludesample.ASVs.excluded.org=M.seq.tax.excludesample.kept.org[rowSums(M.seq.tax.excludesample.kept.org[,1:(ncol(M.seq.tax.excludesample.kept.org)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the original dataset (M.seq.tax))
  M.seq.tax.excludesample.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the subsetted dataset (M.seq.tax.subset))
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.excludesample.excluded.org,M.seq.tax.excludesample.excluded)}
  if (nrow(M.seq.tax.excludesample.ASVs.excluded)!=0&nrow(M.seq.tax.excludesample.ASVs.excluded.org)!=0) { # only happens if at least one ASV was excluded in both dfs
    M.seq.tax.excludesample.kept.org=M.seq.tax.excludesample.kept.org[rowSums(M.seq.tax.excludesample.kept.org[,1:(ncol(M.seq.tax.excludesample.kept.org)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('\nThe following ',nrow(M.seq.tax.excludesample.ASVs.excluded.org),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.excludesample.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.excludesample.kept.org,M.seq.tax.excludesample.ASVs.excluded.org,M.seq.tax.excludesample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.excludesample.ASVs.excluded)!=0&nrow(M.seq.tax.excludesample.ASVs.excluded.org)==0) { # only happens if at least one ASv was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('\nThe following ',nrow(M.seq.tax.excludesample.ASVs.excluded.org),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.excludesample.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.excludesample.kept.org,M.seq.tax.excludesample.ASVs.excluded.org,M.seq.tax.excludesample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.excludesample.ASVs.excluded)==0&nrow(M.seq.tax.excludesample.ASVs.excluded.org)!=0) { # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.excludesample.kept.org=M.seq.tax.excludesample.kept.org[rowSums(M.seq.tax.excludesample.kept.org[,1:(ncol(M.seq.tax.excludesample.kept.org)-ncol.taxo.noNA)])!=0,] # excludes ASVs not longer present in the samples
    cat('\nThe following ',nrow(M.seq.tax.excludesample.ASVs.excluded.org),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(rownames(M.seq.tax.excludesample.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.excludesample.kept.org,M.seq.tax.excludesample.ASVs.excluded.org,M.seq.tax.excludesample.ASVs.excluded)}
  } else if (nrow(M.seq.tax.excludesample.ASVs.excluded)==0&nrow(M.seq.tax.excludesample.ASVs.excluded.org)==0) { # only happens if no ASV was excluded
    cat('\nNo ASVs were excluded due to no observed ASVs in the remaining samples.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.excludesample.kept.org,M.seq.tax.excludesample.ASVs.excluded.org,M.seq.tax.excludesample.ASVs.excluded)}
  }
}  
cat('\n',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.14. Excludes subsetted samples based on the remaining reads ========================================================================================================================

cat('\n\n2.14. Excludes remaining samples based on read counts lower than ',MinimumAllowedReadCount.Analysis,' (MinimumAllowedReadCount.Analysis).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (MinimumAllowedReadCount.Analysis==''|MinimumAllowedReadCount.Analysis==0) {
  cat('\nNo samples were excluded based on their read counts as no MinimumAllowedReadCount.Analysis was chosen or it was set to 0.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
} else {
  M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org=M.seq.tax[,-which(numcolwise(sum)(M.seq.tax)>=MinimumAllowedReadCount.Analysis)] # goes through all numeric columns and excludes those samples (columns) which have more or equal reads (sum) than MinimumAllowedReadCount.Analysis (based on not subsetted df(M.seq.tax))
  M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded=as.data.frame(M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)>=MinimumAllowedReadCount.Analysis)]) # goes through all numeric columns and excludes those samples (columns) which have more or equal reads (sum) than MinimumAllowedReadCount.Analysis (based on  subsetted df(M.seq.tax.subset))
  if (ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)!=ncol.taxo.noNA) { # only if at least one sample was excluded in both dfs
    M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org=M.seq.tax[,-which(numcolwise(sum)(M.seq.tax.subset)<MinimumAllowedReadCount.Analysis)] # only keeps samples with more or equal to MinimumAllowedReadCount reads (based on original df (M.seq.tax))
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)<MinimumAllowedReadCount.Analysis)] # only keeps samples with more or equal to MinimumAllowedReadCount reads (based on subsetted df (M.seq.tax.subset))
    cat('\nThe following ',(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)-ncol.taxo.noNA),' samples were exclulded based on low read counts:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded[,1:((ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded))-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org=M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[rowSums(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[,1:(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the  samples (based on the original dataset (M.seq.tax)
    M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the subsetted dataset (M.seq.tax.subset))
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)}
    if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)!=0&nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)!=0) { # only happens if at least one ASV was excluded in both dfs
      M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org=M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[rowSums(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[,1:(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org)-ncol.taxo.noNA)])!=0,] #excludes ASVs which do not occur anymore in the original df
      M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs which do not occur anymore in the subsetted df
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),' ASVs have been exclulded because they are no longer present in any of the samples:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    } else if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)!=0&nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)==0) { # only happens if at least one ASV was excluded in M.seq.tax.subset but not in M.seq.tax
      M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs which do not occur anymore in the subsetted df
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),' ASVs have been exclulded because they are no longer present in any of the samples:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    } else if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)==0&nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)!=0) { # only happens if at least one ASV was excluded in M.seq.tax but not in M.seq.tax.subset
      M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org=M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[rowSums(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[,1:(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org)-ncol.taxo.noNA)])!=0,] #excludes ASVs which do not occur anymore in the original df
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org),' ASVs have been exclulded because they are no longer present in any of the samples:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    } else if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)==0&nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)==0) { # only happens if no ASV was excluded
      cat('\nNo ASVs were excluded because they do not occur in any of the samples anymore.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    }
  } else if (ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)!=ncol.taxo.noNA&ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)==ncol.taxo.noNA) { # only happens if at least one sample was excluded in the df M.seq.tax.subset but not in M.seq.tax
    M.seq.tax.subset=M.seq.tax.subset[,-which(numcolwise(sum)(M.seq.tax.subset)<MinimumAllowedReadCount.Analysis)] # only keeps samples with more or equal to MinimumAllowedReadCount reads (based on subsetted df (M.seq.tax.subset))
    cat('\nThe following ',(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)-ncol.taxo.noNA),' samples were exclulded based on low read counts:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded[,1:((ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded))-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the subsetted samples (based on the subsetted dataset (M.seq.tax.subset))
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)}
    if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)!=0) { # only happens if at least one ASV was excluded
      M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,] # excludes ASVs which do not occur anymore in the subsetted df
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),' ASVs have been exclulded because they are no longer present in any of the samples:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    } else { # only happens if no ASV was excluded
      cat('\nNo ASVs were excluded because they do not occur in any of the samples anymore.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded)}
    } 
  } else if (ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)==ncol.taxo.noNA&ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)!=ncol.taxo.noNA) { # only happens if at least one sample was excluded in the df M.seq.tax but not in M.seq.tax.subset
    M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org=M.seq.tax[,-which(numcolwise(sum)(M.seq.tax)<MinimumAllowedReadCount.Analysis)] # only keeps samples with more or equal to MinimumAllowedReadCount reads (based on original df (M.seq.tax))
    cat('\nThe following ',(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)-ncol.taxo.noNA),' samples were exclulded based on low read counts:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    cat(colnames(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org[,1:((ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org))-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\n',append=T)
    M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org=M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[rowSums(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[,1:(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org)-ncol.taxo.noNA)])==0,] # creates an ASV by Sample&Taxonomy df with all ASVs which do not occur in the  samples (based on the original dataset (M.seq.tax)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)}
    if (nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)!=0) { # only happens if at least one ASV was excluded
      M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org=M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[rowSums(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org[,1:(ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org)-ncol.taxo.noNA)])!=0,] #excludes ASVs which do not occur anymore in the original df
      cat('The following ',nrow(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org),' ASVs have been exclulded because they are no longer present in any of the samples:\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      cat(rownames(M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=', ',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)}
    } else { # only happens if no ASV was excluded
      cat('\nNo ASVs were excluded because they do not occur in any of the samples anymore.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
      if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.kept.org,M.seq.tax.MinimumAllowedReadCount.Analysis.ASVs.excluded.org)}
    }
  } else if (ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)==ncol.taxo.noNA&ncol(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org)==ncol.taxo.noNA) { # only happens if no sample was excluded in both df.
    cat('\nNo samples were excluded based on low read counts as all samples had more than ',MinimumAllowedReadCount.Analysis,' reads.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
    if (SaveWholeworkspace=='N') {rm(M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded.org,M.seq.tax.MinimumAllowedReadCount.Analysis.samples.excluded)}
  }
}
cat('\n',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# #=== 2.15. Loads and prepares your contextual data ==================================================================================================================================================

cat('\n\n2.15. Load and prepare your contextual data','\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (Metadata =="N"){ # checks whether contextual data was provided and if not it creates a dummy metadata sheet with "A" for Grouping1 and "1" for Grouping2
  M.metadata=as.data.frame(rownames(M.seqtab.nochim))
  row.names(M.metadata)=M.metadata[,1]
  M.metadata[,1]=rep("A",nrow(M.metadata))
  M.metadata$Dummy=rep(1,nrow(M.metadata))
  colnames(M.metadata)[1]="NoMetadata"
  Grouping1="NoMetadata"
  Grouping2="Dummy"
} else if (Metadata=="Y" && Grouping1=="" && Grouping2=="") {
  M.metadata=read.table(file.path(Metadatafilepath),header=T,row.names = NULL)
  rownames(M.metadata)=M.metadata[,1]
  M.metadata=M.metadata[,-1]
  M.metadata$noGroup=rep("A",nrow(M.metadata))
  Grouping1="Group"
  M.metadata$noDummy=rep(1,nrow(M.metadata))
  Grouping2="Dummy"
} else if (Metadata=="Y" && Grouping2=="") {
  M.metadata=read.table(file.path(Metadatafilepath),header=T,row.names = NULL) 
  rownames(M.metadata)=M.metadata[,1]
  M.metadata=M.metadata[,-1]
  M.metadata$Dummy=rep(1,nrow(M.metadata))
  Grouping2="Dummy"
} else {
  M.metadata=read.table(file.path(Metadatafilepath),header=T,row.names = NULL) 
  rownames(M.metadata)=M.metadata[,1]
  M.metadata=M.metadata[,-1]
}

if(is.numeric(M.metadata[,Grouping1])){ # If Grouping1 is a numeric value this adds an A in front of the values.
  M.metadata[,Grouping1]=paste("A",M.metadata[,Grouping1],sep="")
  M.metadata[,Grouping1][M.metadata[,Grouping1]=='ANA']=NA # turns created 'ANA's back to 'NA's
}

cat("The following Paramters are available in your mapfile:\n",file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="",append=TRUE)
cat(colnames(M.metadata),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=", ",append=TRUE)

if (Metadata=='Y'&Grouping1!='NoMetadata') { # checks whether your selected Groupings do occur in your metadata file
  if(Grouping1 %in% colnames(M.metadata)) {
    if (Grouping2 %in% colnames(M.metadata)) {
      cat('\nYour chosen categories (Grouping1: ',Grouping1,', Grouping2: ',Grouping2,') do occur in your metadata file. Yeay!',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="",append=TRUE)
    }
    else {
      cat('\nError: Your Grouping2 does not occur in your Metadata file. Please check spelling.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\t',append=TRUE)
      stop('Your Grouping2 (',Grouping2, ') does not occur in your Metadata file. Please check spelling.')
    }
  } else if (Grouping2 %in% colnames(M.metadata)) {
    cat('\nError: Your Grouping1 (' ,Grouping1,') does not occur in your Metadata file. Please check spelling.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\t',append=TRUE)
    stop("Your Grouping1 (' ,Grouping1,') does not occur in your Metadata file. Please check spelling.")
  } else {
    cat('\nError: None of your Groupings (Grouping1: ' ,Grouping1,', Grouping2: ',Grouping2,') does occur in your Metadata file. Please check spelling.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='\t',append=TRUE)
    stop("None of your Groupings (Grouping1: ' ,Grouping1,', Grouping2: ',Grouping2,') does occur in your Metadata file. Please check spelling.")
  }
}

# #=== 2.16. Excludes Samples based on a category in your map file ==========================================================================================================================

cat('\n\n2.16. Excludes Samples from your metadata sheet based on a category in your map file.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

if (ColumnOfVariableToExclude!='') {
  M.metadata.CategoryToExcluded.samples.excluded=M.metadata[grep(VariableToExclude,M.metadata[,ColumnOfVariableToExclude],invert=F),]   # only takes samples belonging to the VariableToExclude
  M.metadata=M.metadata[grep(VariableToExclude,M.metadata[,ColumnOfVariableToExclude],invert=TRUE),]                                    # Only takes samples not belonging to the VariableToExclude
  if (nrow(M.metadata.CategoryToExcluded.samples.excluded)!=0) {                                                                        # Only happens if at least 1 sample was excluded
    cat('\nThe following ',nrow(M.metadata.CategoryToExcluded.samples.excluded),' samples were excluded because they belong to the Category to exclude (',VariableToExclude,') found in the column ',ColumnOfVariableToExclude,' in your metadata.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat(rownames(M.metadata.CategoryToExcluded.samples.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\n",append=TRUE)
    if (SaveWholeworkspace=='N') {rm(M.metadata.CategoryToExcluded.samples.excluded)}
  } else {
    cat('\nNo samples were excluded because none of the sample in your metadata belong to the Category to exclude (',VariableToExclude,') found in the column ',ColumnOfVariableToExclude,' in your metadata.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    if (SaveWholeworkspace=='N') {rm(M.metadata.CategoryToExcluded.samples.excluded)}
  }
} else {
  cat('\nNo samples will be excluded as no VariableToExclude was chosen. \n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}

# #=== 2.17. Excludes Samples based on selected grouping parameters ===========================================================================================

cat('\n2.17. Excludes Samples based on selected grouping parameters.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
M.metadata.subset=M.metadata
if (SaveWholeworkspace=='N') {rm(M.metadata)}
M.metadata.subset[M.metadata.subset=="NA"]=NA
MapCol1=which(names(M.metadata.subset)==Grouping1) # finds the number of the column of Grouping1 in your metadata
MapCol2=which(names(M.metadata.subset)==Grouping2) # finds the number of the column of Grouping2 in your metadata

cat('\nYour samples will be grouped  based on values in the parameter ',Grouping1,' (Grouping1) and sorted inside these groups based on the parameter ',Grouping2,' (Grouping2).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.metadata.subset.NA=subset(M.metadata.subset,is.na(M.metadata.subset[,MapCol1])|is.na(M.metadata.subset[,MapCol2])) # creates a Sample by Metadata file of samples which have NAs in one of the selected Groupings
if (nrow(M.metadata.subset.NA)==0){ 
  cat('\nNo sample was excluded from your contextual data file due to missing contextual data in the chosen groupings.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
} else {
  cat('\nThe following ',nrow(M.metadata.subset.NA),' samples were excluded from your contextual data sheet because they do not have contextual data in the selected groupings: \n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(rownames(M.metadata.subset.NA),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\n",append=TRUE)
  cat('No metadata in Grouping1: ',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(rownames(subset(M.metadata.subset,is.na(M.metadata.subset[,MapCol1]))),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=", ",append=TRUE)
  cat('\nNo metadata in Grouping2: ',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(rownames(subset(M.metadata.subset,is.na(M.metadata.subset[,MapCol2]))),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=", ",append=TRUE)
}
if (SaveWholeworkspace=='N') {rm(M.metadata.subset.NA)}
M.metadata.subset.noNA=subset(M.metadata.subset,!is.na(M.metadata.subset[,MapCol1])&!is.na(M.metadata.subset[,MapCol2])) # Creates Metadata df without NAs in either of the Grouping columns

# #=== 2.18. Subsets Metadata to fit to your subsetted sequencing data (M.seq.tax.subset) ===================================================================================================

cat('\n\n2.18. Subsets Metadata to fit to your subsetted sequencing data (M.seq.tax.subset).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.metadata.subset=M.metadata.subset.noNA[match((as.character(colnames(M.seq.tax.subset[1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)]))),rownames(M.metadata.subset.noNA)),] # Only keeps samples with the same names as those occuring in the M.seq.tax.subset. This creates empty rows for those samples which do occur in the M.seq.tax.subset but do not occur in M.metadata.subset.noNAs. Does not take the taxonomy in M.seq.tax.subset into account.
M.metadata.subset=subset(M.metadata.subset,!is.na(M.metadata.subset[,MapCol1]),) # deletes newly introduced NA columns

if (nrow(M.metadata.subset.noNA)!=nrow(M.metadata.subset)) { # only happens if at least one sample does not occur in the Metadata which does occur in M.seq.tax.subset and vice versa
  KeepRows=!match(rownames(M.metadata.subset.noNA),(as.character(colnames(M.seq.tax.subset[1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])))) # this is a detour to only keep the samples which do not match to the samples in the M.seq.tax.subset because using !match does result in an empty df
  M.metadata.subset.noseq=cbind(M.metadata.subset.noNA,KeepRows)
  M.metadata.subset.noseq=subset(M.metadata.subset.noseq,is.na(M.metadata.subset.noseq[,ncol(M.metadata.subset.noseq)]))
  M.metadata.subset.noseq=M.metadata.subset.noseq[,-ncol(M.metadata.subset.noseq)]
  cat('\nThe following samples have been excluded from your metadata because they do not occur in your subsetted M.seq.tax.subset: \n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(rownames(M.metadata.subset.noseq),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\n",append=TRUE)
  if (SaveWholeworkspace=='N') {rm(M.metadata.subset.noseq,KeepRows)}
} else {
  cat('\nNo samples have been excluded from your metadata because all samples in your metadata occur in your M.seq.tax.subset.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}

M.metadata.subset.noNA=M.metadata.subset
if (SaveWholeworkspace=='N') {rm(M.metadata.subset)}

# #=== 2.19. Subsets M.seq.tax.subset to samples where metadata is available for the selected grouping parameters. ==========================================================================

cat('\n\n2.19. Subsets M.seq.tax.subset to samples where metadata is available for the selected grouping parameters.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
M.seq.tax.subset.noMetadata=M.seq.tax.subset[,-match(rownames(M.metadata.subset.noNA),(colnames(M.seq.tax.subset)))] # creates table from M.seq.tax.subset only containing samples which do not have metadata available

if(ncol(M.seq.tax.subset.noMetadata)==ncol.taxo.noNA){ # only happens if no samples was excluded (all samples occur in the subsetted metadatafile.)
  cat('\nNo samples have been excluded from your sequencing data because they do not occur in your metadata file (after subsetting).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
} else { # happens if at least one sample was excluded
  M.seq.tax.subset=cbind(M.seq.tax.subset[,match(rownames(M.metadata.subset.noNA),(colnames(M.seq.tax.subset)))],M.seq.tax.subset[,(ncol(M.seq.tax.subset)-ncol.taxo.noNA+1):ncol(M.seq.tax.subset)]) # creates table from M.seq.tax.subset only containing samples which do have metadata available
  cat('\nThe following ',(ncol(M.seq.tax.subset.noMetadata)-ncol.taxo.noNA) ,' samples have been excluded from your sequencing data because they do not occur in your metadatafile or have no data for your selected groupings (after subsetting).\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(colnames(M.seq.tax.subset.noMetadata[,1:(ncol(M.seq.tax.subset.noMetadata)-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\n",append=TRUE)
  cat('\nThe following ',(ncol(M.seq.tax.subset)-ncol.taxo.noNA) ,' samples out of originally ' ,(ncol(M.seq.tax)-ncol.taxo.noNA),'samples remain for further analysis',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.noMetadata)}

M.seq.tax.metadata.ASVs.excluded=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])==0,] # creates ASV by Sample&Taxonomy df with all the ASVs which do not occur anymore in any of the samples

# #=== 2.20 Creates Group and Color vectors to represent the project/samples based on projects/treatments ===================================================================================

cat('\n\n2.20. Creates Group and Color vectors for the visualization.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.projects=as.vector(M.metadata.subset.noNA[,MapCol1]) # vector defining groups/projects based on MapCol1. In original order
M.projects.unique=unique(M.projects) # vector/list of individual projects, in original order
M.projects.unique.ord=sort(M.projects.unique, decreasing=F) # vector/list of individual projects sorted alphabetically

if (M.col=='') {
  M.palette=c("#56B4E9","#E69F00", "#009E73", "#CC79A7","#999999", "#0072B2", "#D55E00","#F0E442") 
  M.col=M.palette[1:length(M.projects.unique)] # This palette is matched to Grouping1 in alphabetical order
  rm(M.palette)
}

M.col.names=(sapply(M.col,color.id)) # gives names for the colors
if (!is.matrix(M.col.names)) {
  M.col.names=ldply (M.col.names, data.frame) # if more than one name for at least one color a list will be created, will be transformed to a data.frame here
  M.col.names=M.col.names[!duplicated((M.col.names$.id)),] # this removes duplicate values
} else {
  M.col.names=t(M.col.names)
}
colnames(M.col.names)=c('HexCode','Name') 

M.match.col=cbind(M.projects.unique.ord,M.col) # table that assigns a specific color to a sample. In alphabetical order
M.match.col.names=cbind(M.projects.unique.ord,M.col.names) # table that assigns a specific color to a sample and shows the name of the color (alphabetical order)

M.groups=match(M.projects, M.match.col) # converts project names to values (in alphabetical order (A=1,B=2....)) This vector is in the original order
M.colvec=mapvalues(M.projects, from = M.match.col[,1], to = M.match.col[,2]) # function in plyr, recodes sample vector. Takes values in M.project and transforms them to the corresponding colors as found in M.match col. This vector is in the original order

M.group.count=matrix(NA,nrow=length(M.projects.unique),ncol=2) # matrix to be filled by the function with as many rows as unique projects and 2 columns
for (i in 1:length(M.projects.unique)) { # Counts how many samples belong to each Group 
  M.group.count[i,1]=M.projects.unique[i]
  M.group.count[i,2]=sum(str_count(M.projects,pattern=M.projects.unique[i])) 
}
M.group.count=as.data.frame(M.group.count) # needs to be transformed to a dataframe to order
M.group.count.ord=M.group.count[order(M.group.count$V1),]  # data frame of individual projects together with number of how often this group appears in the subsorted data sorted alphabetically
M.group.count.ord=as.matrix(M.group.count.ord) # needs to be transformed back to a matrix to work for the txt output.
cat('\nYour chosen groups subsetted for VisuaR contain in total ',length(M.projects),' samples.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
for (i in 1:length(M.projects.unique)) {
  cat('\nGroup ',M.group.count.ord[i,1],' contains ',M.group.count.ord[i,2],' samples and will appear colored in ',M.match.col.names[i,2],'.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
}
if (SaveWholeworkspace=='N') {rm(M.group.count,M.group.count.ord,M.match.col.names)}

# #=== 2.21. Subsampling finished ==========================================================================================================================================================

cat('\n\n2.21. The subsampling is finished.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
if (nrow(M.seq.tax.metadata.ASVs.excluded)!=0) { # Only happens if a sample was excluded and some of the ASVs now do not occur anymore
  M.seq.tax.subset=M.seq.tax.subset[rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)])!=0,]
  cat('\nThe following ',nrow(M.seq.tax.metadata.ASVs.excluded),' ASVs have been excluded because they are no longer present in any of the samples.\n',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat(rownames(M.seq.tax.metadata.ASVs.excluded),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=",",append=TRUE)
} else {
  cat('\nNo ASVs have been excluded because they are no longer present in any of the samples.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}
if (SaveWholeworkspace=='N') {rm(M.seq.tax,M.seq.tax.metadata.ASVs.excluded)}
cat('\nThe following',(ncol(M.seq.tax.subset)-ncol.taxo.noNA),' samples with a total of ',(nrow(M.seq.tax.subset)),' different ASVs out of originally ',nrow(M.seqtab.nochim),' samples with ',(ncol(M.seqtab.nochim)),' different ASVs remain for analysis.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)
cat('\n',colnames(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)]),file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep=',',append=T)
if (SaveWholeworkspace=='N') {rm(M.seqtab.nochim)}

# # Reads of the remaining samples
M.sample.reads.left=colSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)]) # creates a named numeric with the reads per sample
cat('\nThe remaining samples contain in total ',sum(M.sample.reads.left),' ASV counts. With an average of ',round(mean(M.sample.reads.left),1),',  a minimum of ',min(M.sample.reads.left),' and a maximum of ',max(M.sample.reads.left),' ASV counts.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# # ASVs per sample
M.ASV.per.sample=colSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)]!=0) # creates a named numeric with the number of different ASVs per sample
cat('\nThe remaining samples contain in total ',nrow(M.seq.tax.subset),' different ASVs. The samples contain on average ',round(mean(M.ASV.per.sample),1),', with a minimum of ',min(M.ASV.per.sample),' and a maximum of ',max(M.ASV.per.sample),' different ASVs.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

# # Reads of the remaining ASVs
M.ASV.reads.left=rowSums(M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)]) # creates a named numeric with the reads per ASV
cat('\nThe remaining ASVs occur on average ',round(mean(M.ASV.reads.left),1),' times, with a minimum of ',min(M.ASV.reads.left),' times, a maximum of ',max(M.ASV.reads.left),' times.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=T)

closeAllConnections() # closes all currently open connections.

# #=== 3. Community Composition ============================================================================================================================

cat('\n\n\n3. Community Composition',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# #=== 3.1. Creates relative abundance, presence absence tables ============================================================================================================================

cat('\n3.1. Creates presence/absence, relative abundance and other tables.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# # This file includes ASV identities in the first column, followed by respective ASV frequencies, and finally respective ASV taxonomy in the last column
M.seq.tax.subset.print=cbind(rownames(M.seq.tax.subset),data.frame(M.seq.tax.subset,row.names = NULL))
colnames(M.seq.tax.subset.print)[1]="ASV"
write.table(M.seq.tax.subset.print,file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_ASVbySample_abund_final.txt",sep="")),sep='\t',col.names = NA ) #write table
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.print)}

# # Split in two dataframes (Taxonomy and ASVs)
M.seq.tax.subset.ASVs=M.seq.tax.subset[,1:(ncol(M.seq.tax.subset)-ncol.taxo.noNA)] # creates ASV by sample table without taxonomy
M.seq.tax.subset.ASVs=as.data.frame(t(M.seq.tax.subset.ASVs)) 
M.seq.tax.subset.tax=M.seq.tax.subset[,c((ncol(M.seq.tax.subset)-ncol.taxo.noNA+1):ncol(M.seq.tax.subset))] # creates ASV by Taxonomy table

# # merge into one table with relative abundances and Taxonomy
M.ASV.rel=decostand(M.seq.tax.subset.ASVs, method="total") # vegan package, calculates relative abundance
M.seq.tax.subset.rel=cbind(as.data.frame(t(M.ASV.rel)),M.seq.tax.subset.tax) # transforms table and adds taxonomy back to the table
M.seq.tax.subset.rel=cbind(rownames(M.seq.tax.subset.rel),data.frame(M.seq.tax.subset.rel,row.names = NULL)) # modifies table for txt output
colnames(M.seq.tax.subset.rel)[1]="ASV"
write.table(M.seq.tax.subset.rel,file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_ASVbySample_relabund_final.txt",sep="")),sep="\t",col.names = NA)
if (SaveWholeworkspace=='N') {rm(M.ASV.rel,M.seq.tax.subset.rel)}

# # merge into one table with presence absence and taxonomy
M.ASV.pa=decostand(M.seq.tax.subset.ASVs, method="pa") # vegan package, calculates presence absence table
M.seq.tax.subset.pa=cbind(as.data.frame(t(M.ASV.pa)),M.seq.tax.subset.tax) # transforms table and adds taxonomy back to the table
M.seq.tax.subset.pa=cbind(rownames(M.seq.tax.subset.pa),data.frame(M.seq.tax.subset.pa,row.names = NULL)) # modifies table for txt output
colnames(M.seq.tax.subset.pa)[1]="ASV"
write.table(M.seq.tax.subset.pa,file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_ASVbySample_pa_final.txt",sep="")),sep="\t",col.names = NA)
if (SaveWholeworkspace=='N') {rm(M.ASV.pa,M.seq.tax.subset.ASVs,M.seq.tax.subset.pa,M.seq.tax.subset.tax)}

# # Calculates Reads per ASV after subsetting.
M.ASV.reads.left=as.data.frame(M.ASV.reads.left)
M.ASV.reads.left=cbind('ASV'=rownames(M.ASV.reads.left),data.frame(M.ASV.reads.left,row.names = NULL))
colnames(M.ASV.reads.left)[2]="Reads"
write.table(M.ASV.reads.left,file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_ReadsperASV_final.txt",sep="")),sep='\t',col.names = NA )
if (SaveWholeworkspace=='N') {rm(M.ASV.reads.left)}

# #=== 3.2. Creates M.seq.tax.subset.ord: Ordered based on taxonomic levels and based on Grouping ==========================================================================================

# # Orders rows based on taxonomy (A-Z)
ncol.M.seq.tax.subset=ncol(M.seq.tax.subset) # counts how many samples plus taxonomy are in the subsetted M.seq.tax.subset
M.seq.tax.subset.ord=M.seq.tax.subset[order(M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+1)], M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+2)], M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+3)], M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+4)], M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+5)],M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+6)],M.seq.tax.subset[,(ncol.M.seq.tax.subset-ncol.taxo.noNA+7)]),] # Orders the M.seq.tax.subset alphabetically, first based on first taxonomy column, second column.....
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset)}

# #Orders M.seq.tax columns based on projects/groups/contextual data (A-Z or 1-n) 
M.seq.tax.subset.ord.ASV=M.seq.tax.subset.ord[,1:(ncol.M.seq.tax.subset-ncol.taxo.noNA)]                                      # makes ASV table without taxonomy, ordered based on taxonomy
M.seq.tax.subset.ord.tax=M.seq.tax.subset.ord[,c((ncol.M.seq.tax.subset-ncol.taxo.noNA+1):ncol.M.seq.tax.subset)]             # make taxonomy table, alphabetical ordered
M.metadata.subset.noNA.ord=M.metadata.subset.noNA[order(M.metadata.subset.noNA[,MapCol1],M.metadata.subset.noNA[,MapCol2]),]  # orders metadata alphabetically based on the provided Groupings
M.seq.tax.subset.ord.ASV.ord=M.seq.tax.subset.ord.ASV[,rownames(M.metadata.subset.noNA.ord)]                                  # Orders ASV table based on ordered metadata (columns are ordered as they are ordered as rows in M.metadata.subset.noNA.ord)
M.seq.tax.subset.ord=cbind(M.seq.tax.subset.ord.ASV.ord,M.seq.tax.subset.ord.tax)                                             # Rows=ASVs, ordered alphabetically for the taxonomy. Columns=Samples, ordered alphabetically depending on the group they belong to.
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.ord.ASV.ord,M.seq.tax.subset.ord.tax,M.metadata.subset.noNA.ord)}

# #=== 3.3. Creates tables using individual taxonomic levels ===============================================================================================================================

cat('\n\n3.3. Creates tables using individual taxonomic levels',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
CompositionTablePath=file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionTables")
dir.create(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionTables"))
cat('\nThe tables will be saved to ',CompositionTablePath,'.\nFor each taxonomic level 4 tables will be created.\n1. ..TaxonomicLevelBySamples_abund.txt: A Taxonomy by sample table with all reads belonging to the clades summed up and the average of the clades in the last column.\n2. ..TaxonomicLevelBySamples_abund.sort.txt: A Taxonomy by sample table with all reads belonging to the clades summed up and the average of the clades in the last column. Depending on the value set in NumberOfTOPClades the rare taxa will be summed up in the last row and named as others. The Clades are ordered based on descending means.\n3. ..TaxonomicLevelBySamples_relabund_sort.txt: Same as 2. but wiht relative abundances instead of Reads.\n4. ..TaxonomicLevelBySamples_pa.sort.txt: Same as 2. but wiht presences/absences instead of Reads.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# #--- 3.3.1 Phylum level table ----------------------------------------------------------------------------------------------------------------------------------------------------------

# # Phylum level table
M.seq.tax.subset.phylum=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+2))] # Phylum level table; Creates ASV by Sample table only conaining information on Phylum level. This table is ordered by the metadata (alphabetically) and by the taxonomy (alphabetically)

# # Phylum level table with ASVs summed up
M.seq.tax.subset.phylum.sums=aggregate(M.seq.tax.subset.phylum[1:(ncol(M.seq.tax.subset.phylum)-1)], list(Phylum=M.seq.tax.subset.phylum$Phylum),FUN = sum) 

# # Reformat table
row.names(M.seq.tax.subset.phylum.sums)=M.seq.tax.subset.phylum.sums[,1] ;M.seq.tax.subset.phylum.sums=M.seq.tax.subset.phylum.sums[,-1]  # Phylum names are now the rownames, and first column then deleted

# # Add column with row average
M.seq.tax.subset.phylum.sums=cbind(M.seq.tax.subset.phylum.sums,rowMeans(M.seq.tax.subset.phylum.sums)); colnames(M.seq.tax.subset.phylum.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"

# # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
M.seq.tax.subset.phylum.sums.ord=M.seq.tax.subset.phylum.sums[order(-M.seq.tax.subset.phylum.sums$Average),] # Order based on decreasing average
NRP=nrow(M.seq.tax.subset.phylum.sums) # Number of rows in table

# # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
if(NRP>NumberOfTOPClades.box) {
  M.seq.tax.subset.phylum.sums.ord.box=rbind(M.seq.tax.subset.phylum.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.phylum.sums.ord[c((NumberOfTOPClades.box+1):NRP),])); 
  rownames(M.seq.tax.subset.phylum.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" 
} else {
  M.seq.tax.subset.phylum.sums.ord.box=M.seq.tax.subset.phylum.sums.ord
}
# # Calculate relative abundance table
M.seq.tax.subset.phylum.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.phylum.sums.ord.box), method="total"))

# # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
if(NRP>NumberOfTOPClades) {
  M.seq.tax.subset.phylum.sums.ord=rbind(M.seq.tax.subset.phylum.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.phylum.sums.ord[c((NumberOfTOPClades+1):NRP),])); 
  rownames(M.seq.tax.subset.phylum.sums.ord) [(NumberOfTOPClades+1)]="Other" 
}

# # Calculate relative abundances and presence/absence tables
M.seq.tax.subset.phylum.sums.ord.rel=t(decostand(t(M.seq.tax.subset.phylum.sums.ord), method="total"))
M.seq.tax.subset.phylum.sums.ord.pa=t(decostand(t(M.seq.tax.subset.phylum.sums.ord), method="pa"))

write.table(M.seq.tax.subset.phylum.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_PhylumBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.phylum.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_PhylumBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.phylum.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_PhylumBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.phylum.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_PhylumBySamples_abund_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )

if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum.sums,M.seq.tax.subset.phylum.sums.ord,M.seq.tax.subset.phylum.sums.ord.box,M.seq.tax.subset.phylum.sums.ord.pa)}

# #--- 3.3.2. Class level table -----------------------------------------------------------------------------------------------------------------------------------------------------------

# # Class level table
M.seq.tax.subset.class=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+3))]

# # Class level table with ASV summed up
M.seq.tax.subset.class.sums=aggregate(M.seq.tax.subset.class[1:(ncol(M.seq.tax.subset.class)-1)], list(Class=M.seq.tax.subset.class$Class),FUN = sum)

# # Reformat table
row.names(M.seq.tax.subset.class.sums)=M.seq.tax.subset.class.sums[,c(1)]; M.seq.tax.subset.class.sums=M.seq.tax.subset.class.sums[,-c(1)]

# # Add column with row average
M.seq.tax.subset.class.sums=cbind(M.seq.tax.subset.class.sums,rowMeans(M.seq.tax.subset.class.sums)); colnames(M.seq.tax.subset.class.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"

# # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
M.seq.tax.subset.class.sums.ord=M.seq.tax.subset.class.sums[order(-M.seq.tax.subset.class.sums$Average),] # Order based on decreasing average
NRC=nrow(M.seq.tax.subset.class.sums) # Number of rows in table

# # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
if(NRC>NumberOfTOPClades.box) {
  M.seq.tax.subset.class.sums.ord.box=rbind(M.seq.tax.subset.class.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.class.sums.ord[c((NumberOfTOPClades.box+1):NRC),])); 
  rownames(M.seq.tax.subset.class.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
} else {
  M.seq.tax.subset.class.sums.ord.box=M.seq.tax.subset.class.sums.ord
}
# # Calculate relative abundance table
M.seq.tax.subset.class.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.class.sums.ord.box), method="total"))

# # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
if(NRC>NumberOfTOPClades) {
  M.seq.tax.subset.class.sums.ord=rbind(M.seq.tax.subset.class.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.class.sums.ord[c((NumberOfTOPClades+1):NRC),])); 
  rownames(M.seq.tax.subset.class.sums.ord) [(NumberOfTOPClades+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
}

# # Calculate relative abundances and presence/absence tables
M.seq.tax.subset.class.sums.ord.rel=t(decostand(t(M.seq.tax.subset.class.sums.ord), method="total"))
M.seq.tax.subset.class.sums.ord.pa=t(decostand(t(M.seq.tax.subset.class.sums.ord), method="pa"))

write.table(M.seq.tax.subset.class.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_ClassBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.class.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_ClassBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.class.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_ClassBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.class.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_ClassBySamples_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )

if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.class.sums,M.seq.tax.subset.class.sums.ord,M.seq.tax.subset.class.sums.ord.box,M.seq.tax.subset.class.sums.ord.pa)}

# #--- 3.3.3. Order level table -----------------------------------------------------------------------------------------------------------------------------------------------------------

# # Order level table
M.seq.tax.subset.order=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+4))]

# # Order level table with ASV summed up
M.seq.tax.subset.order.sums=aggregate(M.seq.tax.subset.order[1:(ncol(M.seq.tax.subset.order)-1)], list(Order=M.seq.tax.subset.order$Order),FUN = sum)

# # Reformat table
row.names(M.seq.tax.subset.order.sums)=M.seq.tax.subset.order.sums[,c(1)]; M.seq.tax.subset.order.sums=M.seq.tax.subset.order.sums[,-c(1)]

# # Add column with row average
M.seq.tax.subset.order.sums=cbind(M.seq.tax.subset.order.sums,rowMeans(M.seq.tax.subset.order.sums)); colnames(M.seq.tax.subset.order.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"

# # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
M.seq.tax.subset.order.sums.ord=M.seq.tax.subset.order.sums[order(-M.seq.tax.subset.order.sums$Average),] #Order based on decreasing average
NRO=nrow(M.seq.tax.subset.order.sums) # Number of rows in table

# # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
if(NRO>NumberOfTOPClades.box) {
  M.seq.tax.subset.order.sums.ord.box=rbind(M.seq.tax.subset.order.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.order.sums.ord[c((NumberOfTOPClades.box+1):NRO),])); 
  rownames(M.seq.tax.subset.order.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
} else {
  M.seq.tax.subset.order.sums.ord.box=M.seq.tax.subset.order.sums.ord
}
# # Calculate relative abundance table
M.seq.tax.subset.order.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.order.sums.ord.box), method="total"))

# # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
if(NRO>NumberOfTOPClades) {
  M.seq.tax.subset.order.sums.ord=rbind(M.seq.tax.subset.order.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.order.sums.ord[c((NumberOfTOPClades+1):NRO),])); 
  rownames(M.seq.tax.subset.order.sums.ord) [(NumberOfTOPClades+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
}

# # Calculate relative abundances and presence/absence tables
M.seq.tax.subset.order.sums.ord.rel=t(decostand(t(M.seq.tax.subset.order.sums.ord), method="total"))
M.seq.tax.subset.order.sums.ord.pa=t(decostand(t(M.seq.tax.subset.order.sums.ord), method="pa"))

write.table(M.seq.tax.subset.order.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_OrderBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.order.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_OrderBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.order.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_OrderBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.order.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_OrderBySamples_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )

if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.order.sums,M.seq.tax.subset.order.sums.ord,M.seq.tax.subset.order.sums.ord.box,M.seq.tax.subset.order.sums.ord.pa)}

# #--- 3.3.4. Family level table ----------------------------------------------------------------------------------------------------------------------------------------------------------

# # Family level table
M.seq.tax.subset.fam=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+5))]

# # Family level table with ASV summed up
M.seq.tax.subset.fam.sums=aggregate(M.seq.tax.subset.fam[1:(ncol(M.seq.tax.subset.fam)-1)], list(Family=M.seq.tax.subset.fam$Family),FUN = sum)

# # Reformat table
row.names(M.seq.tax.subset.fam.sums)=M.seq.tax.subset.fam.sums[,c(1)]
M.seq.tax.subset.fam.sums=M.seq.tax.subset.fam.sums[,-c(1)]

# # Add column with row average
M.seq.tax.subset.fam.sums=cbind(M.seq.tax.subset.fam.sums,rowMeans(M.seq.tax.subset.fam.sums)); colnames(M.seq.tax.subset.fam.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"

# # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
M.seq.tax.subset.fam.sums.ord=M.seq.tax.subset.fam.sums[order(-M.seq.tax.subset.fam.sums$Average),] #Order based on decreasing average
NRF=nrow(M.seq.tax.subset.fam.sums) # Number of rows in table

# # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
if(NRF>NumberOfTOPClades.box) {
  M.seq.tax.subset.fam.sums.ord.box=rbind(M.seq.tax.subset.fam.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.fam.sums.ord[c((NumberOfTOPClades.box+1):NRF),])); 
  rownames(M.seq.tax.subset.fam.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
} else{
  M.seq.tax.subset.fam.sums.ord.box=M.seq.tax.subset.fam.sums.ord
}
# # Calculate relative abundance table
M.seq.tax.subset.fam.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.fam.sums.ord.box), method="total"))

# # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
if(NRF>NumberOfTOPClades) {
  M.seq.tax.subset.fam.sums.ord=rbind(M.seq.tax.subset.fam.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.fam.sums.ord[c((NumberOfTOPClades+1):NRF),])); 
  rownames(M.seq.tax.subset.fam.sums.ord) [(NumberOfTOPClades+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
}

# # Calculate relative abundances and presence/absence tables
M.seq.tax.subset.fam.sums.ord.rel=t(decostand(t(M.seq.tax.subset.fam.sums.ord), method="total"))
M.seq.tax.subset.fam.sums.ord.pa=t(decostand(t(M.seq.tax.subset.fam.sums.ord), method="pa"))

write.table(M.seq.tax.subset.fam.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_FamilyBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.fam.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_FamilyBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.fam.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_FamilyBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.fam.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_FamilyBySamples_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )

if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.fam.sums,M.seq.tax.subset.fam.sums.ord,M.seq.tax.subset.fam.sums.ord.box,M.seq.tax.subset.fam.sums.ord.pa)}

# #--- 3.3.5. Genus level table -----------------------------------------------------------------------------------------------------------------------------------------------------------

# # Genus level table
M.seq.tax.subset.genus=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+6))]

# # Genus level table with ASV summed up
M.seq.tax.subset.genus.sums=aggregate(M.seq.tax.subset.genus[1:(ncol(M.seq.tax.subset.genus)-1)], list(Genus=M.seq.tax.subset.genus$Genus),FUN = sum)

# # Reformat table
row.names(M.seq.tax.subset.genus.sums)=M.seq.tax.subset.genus.sums[,c(1)]; M.seq.tax.subset.genus.sums=M.seq.tax.subset.genus.sums[,-c(1)]

# # Add column with row average
M.seq.tax.subset.genus.sums=cbind(M.seq.tax.subset.genus.sums,rowMeans(M.seq.tax.subset.genus.sums)); colnames(M.seq.tax.subset.genus.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"

# # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
M.seq.tax.subset.genus.sums.ord=M.seq.tax.subset.genus.sums[order(-M.seq.tax.subset.genus.sums$Average),] #Order based on decreasing average
NRG=nrow(M.seq.tax.subset.genus.sums) # Number of rows in table

# # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
if(NRG>NumberOfTOPClades.box) {
  M.seq.tax.subset.genus.sums.ord.box=rbind(M.seq.tax.subset.genus.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.genus.sums.ord[c((NumberOfTOPClades.box+1):NRG),])); 
  rownames(M.seq.tax.subset.genus.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
} else {
  M.seq.tax.subset.genus.sums.ord.box=M.seq.tax.subset.genus.sums.ord
}
# # Calculate relative abundance table
M.seq.tax.subset.genus.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.genus.sums.ord.box), method="total"))

# # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
if(NRG>NumberOfTOPClades) {
  M.seq.tax.subset.genus.sums.ord=rbind(M.seq.tax.subset.genus.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.genus.sums.ord[c((NumberOfTOPClades+1):NRG),])); 
  rownames(M.seq.tax.subset.genus.sums.ord) [(NumberOfTOPClades+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
}

# # Calculate relative abundances and presence/absence tables
M.seq.tax.subset.genus.sums.ord.rel=t(decostand(t(M.seq.tax.subset.genus.sums.ord), method="total"))
M.seq.tax.subset.genus.sums.ord.pa=t(decostand(t(M.seq.tax.subset.genus.sums.ord), method="pa"))

write.table(M.seq.tax.subset.genus.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_GenusBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.genus.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_GenusBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.genus.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_GenusBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
write.table(M.seq.tax.subset.genus.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_GenusBySamples_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )

if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.genus.sums,M.seq.tax.subset.genus.sums.ord,M.seq.tax.subset.genus.sums.ord.box,M.seq.tax.subset.genus.sums.ord.pa)}

# #--- 3.3.6. Species level table ---------------------------------------------------------------------------------------------------------------------------------------------------------

if (ncol.taxo.noNA==7) {
  # # Species Level table
  M.seq.tax.subset.species=M.seq.tax.subset.ord[,c(1:(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA),(ncol(M.seq.tax.subset.ord)-ncol.taxo.noNA+7))]
  
  # # Species level table with ASV summed up
  M.seq.tax.subset.species.sums=aggregate(M.seq.tax.subset.species[1:(ncol(M.seq.tax.subset.species)-1)], list(Species=M.seq.tax.subset.species$Species),FUN = sum)
  
  # # Reformat table
  row.names(M.seq.tax.subset.species.sums)=M.seq.tax.subset.species.sums[,1]; M.seq.tax.subset.species.sums=M.seq.tax.subset.species.sums[,-1]
  
  # # Add column with row average
  M.seq.tax.subset.species.sums=cbind(M.seq.tax.subset.species.sums,rowMeans(M.seq.tax.subset.species.sums)); colnames(M.seq.tax.subset.species.sums) [ncol.M.seq.tax.subset-ncol.taxo.noNA+1]="Average"
  
  # # Sort table based on decreasing average, show top 20 taxa and combine other rare taxa to "Other"
  M.seq.tax.subset.species.sums.ord=M.seq.tax.subset.species.sums[order(-M.seq.tax.subset.species.sums$Average),] #Order based on decreasing average
  NRS=nrow(M.seq.tax.subset.species.sums) # Number of rows in table
  
  # # Takes the first 'NumberOfTOPClades.box' rows; sums other rows with rare taxa
  if(NRS>NumberOfTOPClades.box) {
    M.seq.tax.subset.species.sums.ord.box=rbind(M.seq.tax.subset.species.sums.ord[c(1:NumberOfTOPClades.box),],colSums(M.seq.tax.subset.species.sums.ord[c((NumberOfTOPClades.box+1):NRG),])); 
    rownames(M.seq.tax.subset.species.sums.ord.box) [(NumberOfTOPClades.box+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
  } else {
    M.seq.tax.subset.species.sums.ord.box=M.seq.tax.subset.species.sums.ord
  }
  # # Calculate relative abundance table
  M.seq.tax.subset.species.sums.ord.box.rel=t(decostand(t(M.seq.tax.subset.species.sums.ord.box), method="total"))
  
  # # Takes the first 'NumberOfTOPClades' rows; sums other rows with rare taxa
  if(NRS>NumberOfTOPClades) {
    M.seq.tax.subset.species.sums.ord=rbind(M.seq.tax.subset.species.sums.ord[c(1:NumberOfTOPClades),],colSums(M.seq.tax.subset.species.sums.ord[c((NumberOfTOPClades+1):NRG),])); 
    rownames(M.seq.tax.subset.species.sums.ord) [(NumberOfTOPClades+1)]="Other" #take first 20 rows; sum other rows with rare taxa; call new row "Other"
  }
  
  #Calculate relative abundances and presence/absence tables
  M.seq.tax.subset.species.sums.ord.rel=t(decostand(t(M.seq.tax.subset.species.sums.ord), method="total"))
  M.seq.tax.subset.species.sums.ord.pa=t(decostand(t(M.seq.tax.subset.species.sums.ord), method="pa"))
  
  write.table(M.seq.tax.subset.species.sums,file.path(CompositionTablePath,paste(VisuaRProjectName,"_SpeciesBySamples_abund.txt",sep="")),sep='\t',col.names = NA )
  write.table(M.seq.tax.subset.species.sums.ord,file.path(CompositionTablePath,paste(VisuaRProjectName,"_SpeciesBySamples_abund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
  write.table(M.seq.tax.subset.species.sums.ord.rel,file.path(CompositionTablePath,paste(VisuaRProjectName,"_SpeciesBySamples_relabund_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
  write.table(M.seq.tax.subset.species.sums.ord.pa,file.path(CompositionTablePath,paste(VisuaRProjectName,"_SpeciesBySamples_pa_TOP",NumberOfTOPClades,".txt",sep="")),sep='\t',col.names = NA )
  
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.species.sums,M.seq.tax.subset.species.sums.ord,M.seq.tax.subset.species.sums.ord.box,M.seq.tax.subset.species.sums.ord.pa)}
}
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.ord)}

# #=== 3.4. Creates Bubble Plots and BoxPlots of relative sequence abundances for all taxonomic levels ==============================================================================================

CompositionPlotsPath=file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots") 
dir.create(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots"))          # Creates Folder for composition plots

cat('\n\n3.4. Create Bubble and BoxPlots of relative abundances for all taxonomic levels',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nThe plots will be saved to ',CompositionPlotsPath,'.\nFor each taxonomic level a seperated folder will be created inside this folder.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\n1. The bubbleplot will show the samples on the x-axis (grouped and alphabetically sorted by ',Grouping1,' (Grouping1) and inside the groups sorted on increasing ',Grouping2,' (Grouping2)). The TOP',NumberOfTOPClades,' clades will be shown on the y axis. Relative Abundances will be presented as dot size.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\n2. The boxplot BoxByClade_RA will show the TOP',NumberOfTOPClades.box,' clades on the x-axis and the relative abundance on the y axis. For each Clade ',length(M.projects.unique),' boxplots will be shown (Number of your Groups provided as Grouping1).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\n3. The boxplot BoxByGroup_RA will show the ',length(M.projects.unique),' Groups provided as Grouping1 on the x axis and the relative abundance on the y axis. For each Clade 1 boxplot will be shown per Group.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\n4. ',NumberOfTOPClades,' Boxplots will be produced of the TOP',NumberOfTOPClades ,' clades. Clade_1 is the Clade with the highest abundance, Clade_2 the one with the second highest abundance and so on.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

TVP=M.seq.tax.subset.phylum.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.phylum.sums.ord.rel)-1))] # takes a ASV x Sample Table. In this table the samples are ordered by the Groupings (1. Grouping1, 2. Grouping2)
TVC=M.seq.tax.subset.class.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.class.sums.ord.rel)-1))]
TVO=M.seq.tax.subset.order.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.order.sums.ord.rel)-1))]
TVF=M.seq.tax.subset.fam.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.fam.sums.ord.rel)-1))]
TVG=M.seq.tax.subset.genus.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.genus.sums.ord.rel)-1))]

TVP.b=M.seq.tax.subset.phylum.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.phylum.sums.ord.box.rel)-1))] #takes a ASV x Sample Table, removes last column (average)
TVC.b=M.seq.tax.subset.class.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.class.sums.ord.box.rel)-1))] 
TVO.b=M.seq.tax.subset.order.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.order.sums.ord.box.rel)-1))] 
TVF.b=M.seq.tax.subset.fam.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.fam.sums.ord.box.rel)-1))] 
TVG.b=M.seq.tax.subset.genus.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.genus.sums.ord.box.rel)-1))] 


if (ncol.taxo.noNA==7) {
  TVS=M.seq.tax.subset.species.sums.ord.rel[,c(1:(ncol(M.seq.tax.subset.species.sums.ord.rel)-1))]
  TVS.b=M.seq.tax.subset.species.sums.ord.box.rel[,c(1:(ncol(M.seq.tax.subset.species.sums.ord.box.rel)-1))]
  levels=list(TVP,TVC,TVO,TVF,TVG,TVS)
  levels.b=list(TVP.b,TVC.b,TVO.b,TVF.b,TVG.b,TVS.b)
  tax.levels=c("Phylum","Class","Order","Family","Genus","Species")
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum.sums.ord.rel,M.seq.tax.subset.class.sums.ord.rel,M.seq.tax.subset.order.sums.ord.rel,M.seq.tax.subset.fam.sums.ord.rel,M.seq.tax.subset.genus.sums.ord.rel,M.seq.tax.subset.species.sums.ord.rel)}
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum.sums.ord.box.rel,M.seq.tax.subset.class.sums.ord.box.rel,M.seq.tax.subset.order.sums.ord.box.rel,M.seq.tax.subset.fam.sums.ord.box.rel,M.seq.tax.subset.genus.sums.ord.box.rel,M.seq.tax.subset.species.sums.ord.box.rel)}
} else {
  levels=list(TVP,TVC,TVO,TVF,TVG)
  levels.b=list(TVP.b,TVC.b,TVO.b,TVF.b,TVG.b)
  tax.levels=c("Phylum","Class","Order","Family","Genus")
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum.sums.ord.rel,M.seq.tax.subset.class.sums.ord.rel,M.seq.tax.subset.order.sums.ord.rel,M.seq.tax.subset.fam.sums.ord.rel,M.seq.tax.subset.genus.sums.ord.rel)}
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum.sums.ord.box.rel,M.seq.tax.subset.class.sums.ord.box.rel,M.seq.tax.subset.order.sums.ord.box.rel,M.seq.tax.subset.fam.sums.ord.box.rel,M.seq.tax.subset.genus.sums.ord.box.rel)}
}

M.projects.ord=sort(M.projects) # sorts the color vector alphabetically to fit TVP, TVC etc.
M.colvec.ord=mapvalues(M.projects.ord, from = M.match.col[,1], to = M.match.col[,2]) # function in library (plyr), recodes sample vector, creates a color vector for the alphabetically by Grouping sorted data

# # only taxonomic levels with more than 1 lineage will be visualized. 
if (is.matrix(TVS)&is.matrix(TVS.b)) {
  startclade=6
  if (is.matrix(TVG)&is.matrix(TVG.b)) {
    startclade=5
    if (is.matrix(TVF)&is.matrix(TVF.b)) {
      startclade=4
      if (is.matrix(TVO)&is.matrix(TVO.b)) {
        startclade=3
        if (is.matrix(TVC)&is.matrix(TVC.b)) {
          startclade=2
          if (is.matrix(TVP)&is.matrix(TVP.b)) {
            startclade=1
            cat('\nPlots will be produced for all taxonomic levels.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
          } else {
            cat('\nNo plots will be produced on phylum level because only one phylum is present after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
          }
        } else {
          cat('\nNo plots will be produced on phylum level because only one phylum is present after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
          cat('\nNo plots will be produced on class level because only one class is present after subsetting (',as.character(M.seq.tax.subset.class[1,ncol(M.seq.tax.subset.class)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
        }
      } else {
        cat('\nNo plots will be produced on phylum level because only one phylum is present after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
        cat('\nNo plots will be produced on class level because only one class is present after subsetting (',as.character(M.seq.tax.subset.class[1,ncol(M.seq.tax.subset.class)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
        cat('\nNo plots will be produced on order level because only one order is present after subsetting (',as.character(M.seq.tax.subset.order[1,ncol(M.seq.tax.subset.order)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
      }
    } else {
      cat('\nNo plots will be produced on phylum level because only one phylum is present after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
      cat('\nNo plots will be produced on class level because only one class is present after subsetting (',as.character(M.seq.tax.subset.class[1,ncol(M.seq.tax.subset.class)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
      cat('\nNo plots will be produced on order level because only one order is present after subsetting (',as.character(M.seq.tax.subset.order[1,ncol(M.seq.tax.subset.order)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
      cat('\nNo plots will be produced on family level because only one family is present after subsetting (',as.character(M.seq.tax.subset.fam[1,ncol(M.seq.tax.subset.fam)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    }
  } else {
    cat('\nNo plots will be produced on phylum level because only one phylum is present after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat('\nNo plots will be produced on class level because only one class is present after subsetting (',as.character(M.seq.tax.subset.class[1,ncol(M.seq.tax.subset.class)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat('\nNo plots will be produced on order level because only one order is present after subsetting (',as.character(M.seq.tax.subset.order[1,ncol(M.seq.tax.subset.order)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat('\nNo plots will be produced on family level because only one family is present after subsetting (',as.character(M.seq.tax.subset.fam[1,ncol(M.seq.tax.subset.fam)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
    cat('\nNo plots will be produced on genus level because only one genus is present after subsetting (',as.character(M.seq.tax.subset.genus[1,ncol(M.seq.tax.subset.genus)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  }
} else {
  startclade=7
  cat('\nNo plots will be produced because all levels contain only one clade after subsetting (',as.character(M.seq.tax.subset.phylum[1,ncol(M.seq.tax.subset.phylum)]),').',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}

if (ncol.taxo.noNA==7) {
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum,M.seq.tax.subset.class,M.seq.tax.subset.order,M.seq.tax.subset.fam,M.seq.tax.subset.genus,M.seq.tax.subset.species)}
} else {
  if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.phylum,M.seq.tax.subset.class,M.seq.tax.subset.order,M.seq.tax.subset.fam,M.seq.tax.subset.genus)}
}

Sys.sleep(10) # this stops the execution for 10 seconds. If this is not set, the for loop is not always correctly executed

for (k in startclade:(ncol.taxo.noNA-1)){ 
  # # Reads the matrix (1 for the bubble plots and 1 for the box plots)
  TV=as.matrix(as.data.frame(levels[k]))      # takes the respective relative abundance matrix for the bubble plots and the single-clade box plots (Clade by Sample)
  TV.b=as.matrix(as.data.frame(levels.b[k]))  # takes the respective relative abundance matrix for the grouped box plots (Clade by Sample)
  nrow.TV=nrow(TV)                            # number of rows = number of different clades
  nrow.TV.b=nrow(TV.b)                        # number of rows = number of different clades
  
  # # Creates group file for coloring of the plot.
  TV.col.group=c(rep(M.projects,nrow.TV))       # group file: fits to concatenated sample column, because the order of Groups (colnames) is multiplied by nrow.TV
  TV.col.group.b=c(rep(M.projects,nrow.TV.b))
  TV.col.group=sort(TV.col.group, decreasing=F) # sorts new group file alphabetically/increasing
  TV.col.group.b=sort(TV.col.group.b,decreasing = F)
  
  TV[TV==0]=NA # turns zeros in relative abundance file into NAs, elsewise zeros are shown as dots in the bubble plot
  TV.b[TV.b==0]=NA
  
  # # prepares relative sequence abundance tables for ggplot
  TV.levels=rownames(TV)                                                          # Level IDs (clade names) used for bubbleplot
  TVX <- melt(TV, id.vars = "Clade", variable.name="Sample", value.name = "Size") # turns multiple "Sample" columns (x axis) into one concatenated "Sample" column 
  TVX.b=melt(TV.b,id.vars="Clade", variable.name="Sample", value.name = "Size")
  colnames(TVX)=c("Clade","Sample","Size")
  colnames(TVX.b)=c("Clade","Sample","Size")
  
  dir.create(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots",tax.levels[k])) 
  
  # #---Create bubble plots 
  pdf(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots",tax.levels[k],paste(VisuaRProjectName,"_Bubble_RA.pdf",sep="")),height=(NumberOfTOPClades*0.3),width=((ncol(TV)*0.2)+((k+1)*1.7)),useDingbats=F) # useDingbats=F is needed to avoid font display problems with downstream graphics software
  M.bubble=ggplot(TVX, aes(x = Sample, y = factor(TVX$Clade, levels=rev(TV.levels)), col=factor(TV.col.group))) +   # plots samples vs Clades using colors according to group file (TV.col group is ordered alphabetically/increasing))
    geom_point(aes(size = Size*100)) +                                                                              # adjust size of circles for visualization reasons
    scale_colour_manual(values=M.col) +                                                                             # the same colors for projects throughout workflow. M.col has as many colors as groups and the first color in M.col is assigned to the first group in alphabetical/increasing order
    labs(x="Samples", y="Clades", col=Grouping1,size="Relative\nSequence\nAbundance\n(in %)",title=paste("Relative Sequence Abundance, ",tax.levels[k]," level ","(",VisuaRProjectName,")",sep="")) + #axis labels
    theme_bw() +
    theme(axis.text.x = element_text(face="plain", color="Black", size=10, angle=90,hjust=1), 
          axis.text.y = element_text(face="plain", color="Black", size=10, angle=0),
          plot.title = element_text(size=8))
  print(M.bubble) # the plot needs to be printed if the script is run at once. Elsewise an empty plot will be created.
  dev.off()
  
  # #---Creates boxplot of top clades. Main Group=Clade, SubGroup=Grouping1
  TVX.box=cbind(TVX.b,TV.col.group.b) # attaches the group column to the concatenated table
  TVX.box$Size[is.na(TVX.box$Size)]=0
  pdf(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots",tax.levels[k],paste(VisuaRProjectName,"_",tax.levels[k],"_TOP",NumberOfTOPClades.box,"_BoxByClade_RA.pdf",sep="")),height=(4+(k*2)),width=NumberOfTOPClades.box*1.2,useDingbats=F) 
  M.box=ggplot(TVX.box,aes(x=Clade,y=Size*100,fill=factor(TV.col.group.b),colour="black"))+ # Transforms relative sequence Abundance to % values
    stat_boxplot(geom = "errorbar",position =position_dodge(preserve='single'))+            
    geom_boxplot(varwidth=F, notch=F, outlier.size=1,outlier.colour = "black",position=position_dodge(preserve = "single")) + 
    labs(x="Clades",y=paste("Relative Sequence Abundance [%]"),title=paste("Relative Sequence Abundance, ",tax.levels[k]," level ","(",VisuaRProjectName,")",sep="")) +
    scale_color_manual(values = "black",guide=FALSE) +
    scale_fill_manual(drop=F,values=M.col) + 
    theme(axis.text = element_text(colour="black"),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x=element_text(angle=90,colour="black",size=14), #,vjust=-0.1
          axis.text.y=element_text(colour="black",size=14),
          axis.title = element_text(colour="black",size=16),
          axis.line=element_line(colour="black"),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(colour="grey",size=0.5),
          panel.grid.minor.y = element_line(colour="grey",size = 0.25),
          panel.border = element_rect(fill=NA,colour = "black",size=0.5),
          legend.box.background = element_rect(colour="black"),
          legend.box.margin = margin(0.18,0.1,0.1,0.5),
          legend.key=element_rect(fill="white"),
          legend.text=element_text(size=14),
          axis.ticks=element_line(colour = "black"))
  print(M.box)
  dev.off()
  
  # #---Create boxplots of top clades Main Group=Grouping1, SubGroup=Clade, Colors:Rainbow
  pdf(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots",tax.levels[k],paste(VisuaRProjectName,"_BoxByGroup_RA.pdf",sep="")),height=(6),width=(NumberOfTOPClades*0.5)+((k-1)*2.5),useDingbats=F) #height=(8),width=(25)
  M.box.2=ggplot(TVX.box,aes(x=TV.col.group.b,y=Size*100,fill=(Clade),colour='black')) +
    stat_boxplot(geom = "errorbar",position =position_dodge(preserve='single')) + 
    geom_boxplot(varwidth=F, notch=F, outlier.size=1,position=position_dodge(preserve = "single")) + 
    labs(x=Grouping1,y=paste("Relative Sequence Abundance [%]"),title=paste("Relative Sequence Abundance, ",tax.levels[k]," level ","(",VisuaRProjectName,")",sep="")) +
    scale_x_discrete()+ # Orders the groups in the plot
    scale_color_manual(values = "black",guide=FALSE)+
    scale_fill_brewer(palette = 'Paired') + 
    guides(fill=guide_legend(ncol=1))+ 
    theme(axis.text = element_text(colour="black"),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x=element_text(angle=0,colour="black",size=14),
          axis.text.y=element_text(colour="black",size=14),
          axis.title = element_text(colour="black",size=16),
          axis.line=element_line(colour="black"),
          panel.background = element_blank(),
          panel.grid.major.y = element_line(colour="grey",size=0.5),
          panel.grid.minor.y = element_line(colour="grey",size = 0.25),
          panel.border = element_rect(fill=NA,colour = "black",size=0.5),
          legend.box.background = element_rect(colour="black"),
          legend.box.margin = margin(0.18,0.1,0.1,0.5),
          legend.key=element_rect(fill="white"),
          legend.text=element_text(size=14),
          axis.ticks=element_line(colour = "black"))
  print(M.box.2)
  dev.off()
  
  # #---Makes boxplots of all Top Clades seperately
  if (length(M.projects.unique)>1){               # only if more than 1 group in Grouping1
    BoxClades=unique(TVX[,1])                     # list of unique clades in TVX
    for(f in seq_along(BoxClades)) {              # for-loop goes through every single clade
      BoxClade=as.character(BoxClades[f])         # this clade will be plotted
      TVX.single.box=TVX[grep(BoxClade,TVX[,1]),] # takes all relative abundance values belonging to the selected clade, adds a column with Grouping1, replaces NAs with 0
      TVX.single.box.groups=cbind(TVX.single.box,M.projects.ord)
      TVX.single.box.groups[is.na(TVX.single.box.groups)]=0
      pdf(file.path(PathToVisuaRAnalysis,"Alpha_Diversity","CompositionPlots",tax.levels[k],paste("Clade_",f,".pdf",sep="")),height = 4,width=length(M.projects.unique),useDingbats=F) # height=3,width=((2)*length(M.projects.unique))
      M.single.box=ggplot(TVX.single.box.groups, aes(x=M.projects.ord,y=Size*100,group=M.projects.ord))
      M.single.box= M.single.box +
        stat_boxplot(geom = "errorbar", width = 0.1) +                  
        geom_boxplot(varwidth=F, notch=F, fill=M.col,outlier.size=1) +  
        geom_signif(comparisons=combn(sort(as.character(unique(TVX.single.box.groups$M.projects.ord))),2,simplify=F),
                    step_increase=0.1,
                    size=0.3,
                    textsize=3,
                    tip_length=0.01,
                    map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))+
        stat_n_text(size=4)+
        labs(y="Relative Sequence Abundance [%]",title=BoxClade) +
        stat_summary(fun.y=mean,colour="black",geom="point",shape=17) +
        theme(axis.text = element_text(colour="black"),
              legend.title = element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.x=element_text(colour="black",size=10),
              axis.text.y=element_text(colour="black",size=12),
              # axis.title.x = element_text(colour="black",size=14),
              axis.title.x=element_blank(),
              axis.title.y = element_text(colour='black',size=12),
              plot.title = element_text(colour='black',size = 3),
              axis.line=element_line(colour="black"),
              panel.background = element_blank(),
              panel.grid.major.y = element_line(colour="grey",size=0.25),
              panel.grid.minor.y = element_blank(),
              panel.border = element_rect(fill=NA,colour = "black",size=0.5),
              legend.box.background = element_rect(colour="black"),
              legend.box.margin = margin(0.18,0.1,0.1,0.5),
              legend.key=element_rect(fill="white"),
              legend.text=element_text(size=12),
              axis.ticks=element_line(colour = "black"))
      print(M.single.box)
      dev.off()
    }
  }
}

if (ncol.taxo.noNA==7) {
  if (SaveWholeworkspace=='N') {rm(TV,TV.b,levels,levels.b,M.box,M.box.2,M.bubble,M.single.box,TVC,TVC.b,TVF,TVF.b,TVG,TVG.b,TVO,TVO.b,TVP,TVP.b,TVS,TVS.b,TVX,TVX.b,TVX.box,TVX.single.box,TVX.single.box.groups)}
} else {
  if (SaveWholeworkspace=='N') {rm(TV,TV.b,levels,levels.b,M.box,M.box.2,M.bubble,M.single.box,TVC,TVC.b,TVF,TVF.b,TVG,TVG.b,TVO,TVO.b,TVP,TVP.b,TVX,TVX.b,TVX.box,TVX.single.box,TVX.single.box.groups)}
}

closeAllConnections() # closes all currently open connections.

# #================== 4. Alpha Diversity ================================================================================================================================================

cat('\n\n4. Alpha Diversity',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M=t(M.seq.tax.subset.ord.ASV)       # Sample by ASV table ordered by taxonomy. Not ordered by Grouping1.
if (SaveWholeworkspace=='N') {rm(M.seq.tax.subset.ord.ASV)}
# write.table(M,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_ASVbySample_abund.txt",sep="")),sep='\t',col.names = NA) 
M.rel=decostand(M, method="total")  # Calculates relative sequence abundances
M.pa=decostand(M, method="pa")      # Calculates presence/absence

M.sobs=rowSums(M.pa)  # Number of different observed ASVs per sample
M.ASVs=rowSums(M)     # Number of reads per sample, not ordered based on Grouping1

# #=======  4.1. Calculates diversity indices ============================================================================================================================================

cat('\n4.1. Calculates diversity indices',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.minreads=min(M.ASVs) # minimum number of reads in any given sample. Will be used to subsample in alpha diversity indices calculation. This allows comparison of samples of different sequencing depth.
cat('\nThe alpha diversity indices will be calculated with a subsample of ',M.minreads,' ASVs (minimum number of reads occuring in any sample).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nObserved ASV richness, Chao1 richness estimate, Shannon Entropy, Inverse Simpson Diversity, absolute and relative Singletons can be found in the file ',file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_DiversityIndices_",Grouping1,".txt",sep="")),'.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nNote: diversity indices may be based on iterative subsampling, whereas SSOabs and SSOrel are based on the whole dataset.\nThis can result in the odd case that there seem to be more observed relative ASVs than observed ASVs.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

NR=nrow(M) # number of rows of matrix M = number of Samples

# #--- 4.1.1. Observed ASV richness subsampled -------------------------------------------------------------------------------------------------------------------------------------------

M.sobs.r=matrix(NA,nrow=NR,ncol=NI)                                                   # matrix to be filled by the function

for (j in 1:NR) {                                                                     # calculates iterated observed ASV richness, subsampled using M.minreads
  for (i in 1:NI) {
    M.sobs.r[j,i]=estimateR(rrarefy(M[j,],sample=M.minreads),replace=F)[1,1]          # replace=T, if sample-parameter is larger than the smallest sample (not recommended)
  } 
}

M.sobs.r.mean=rowMeans(M.sobs.r)                                                      # calculation of average observed ASV richness from (NI reps)

if (SaveWholeworkspace=='N') {rm(M.sobs.r)}

# #--- 4.1.2. Chao1 Richness Estimate ----------------------------------------------------------------------------------------------------------------------------------------------------

M.chao1.r=matrix(NA, nrow=NR,ncol=NI)                                                 # matrix to be filled by the function

for (j in 1:NR) {                                                                     # calculates iterated Chao1 richness estimate, subsampled using M.minreads
  for (i in 1:NI) {
    M.chao1.r[j,i]=estimateR(rrarefy(M[j,],sample=M.minreads),replace=F)[2,1] 
  }
}

M.chao1.r.mean=rowMeans(M.chao1.r)                                                    #calculation of average Chao richness estimate from (NI reps)

# #--- 4.1.3. Shannon Entropy ------------------------------------------------------------------------------------------------------------------------------------------------------------

M.shan.r=matrix(NA,nrow=NR,ncol=NI)                                                   # matrix to be filled by the function

for (j in 1:NR) {                                                                     # calculate iterated Shannon entropy index, subsampled using M.minreads
  for (i in 1:NI) {                                                                   # happens once for each iteration
    M.shan.r[j,i]=diversity(rrarefy(M[j,],sample=M.minreads), index="shannon") 
  }
}      

M.shan.r.mean=rowMeans(M.shan.r)                                                      # calculation of average Shannon Entropy from (NI reps)

if (SaveWholeworkspace=='N') {rm(M.shan.r)}

# #--- 4.1.4. Inverse Simpson Diversity --------------------------------------------------------------------------------------------------------------------------------------------------

M.invs.r=matrix(NA,nrow=NR,ncol=NI)                                                   # matrix to be filled by the function

for (j in 1:NR) {                                                                     # calculates iterated inverse simpson diversity index, subsampled using M.minreads
  RE=rep(NA,NI) 
  for (i in 1:NI) {
    RE[i]=diversity(rrarefy(M[j,],sample=M.minreads), index="invsimpson") 
  }
  M.invs.r[j,]= RE;
}     

M.invs.r.mean=rowMeans(M.invs.r)                                                      # calculation of average Inverse Simpson Diversity from (NI reps)

if (SaveWholeworkspace=='N') {rm(M.invs.r)}

# #--- 4.1.5. Number of Absolute Singletons ----------------------------------------------------------------------------------------------------------------------------------------------

M.SSASV.abs=as.matrix(M[,apply(M,2,sum)==1])                                          # Table of only absolute singletons, 2 for columns

if (is.null(ncol(M.SSASV.abs))){                                                      # only happens if data set has no absolute Single Sequence ASV (SSASVabs), i.e. absolute singletons
  M.SSASV.abs.nr=rep(0,NR); names(M.SSASV.abs.nr)=names(M.sobs)                       # creates a string of 0's (length=NR=number of samples)
  M.SSASV.abs.pc=rep(0,NR); names(M.SSASV.abs.pc)=names(M.sobs)
  M.SSASV.abs.nr.t=t(M.SSASV.abs.nr); row.names(M.SSASV.abs.nr.t)='M.SSASV.abs.nr'    # transposed to fit combined table below
  M.SSASV.abs.pc.t=t(M.SSASV.abs.pc); row.names(M.SSASV.abs.pc.t)="M.SSASV.abs.pc"    # transposed to fit combined table below
  if (SaveWholeworkspace=='N') {rm(M.SSASV.abs,M.SSASV.abs.nr,M.SSASV.abs.pc)}
} else if (ncol(M.SSASV.abs)==1) {
  M.SSASV.abs.name=which(colSums(M)==1)                                               # finds name and position of absolute singleton in M
  colnames(M.SSASV.abs)=as.character(names(M.SSASV.abs.name))                         # renames the ASV to actual ASV name
  M.SSASV.abs.nr=rowSums(M.SSASV.abs)                                                 # Number of SSASVabs per sample
  M.SSASV.abs.pc=round(M.SSASV.abs.nr/M.sobs*100,1)                                   # % SSASVabs
  M.SSASV.abs.nr.t=t(M.SSASV.abs.nr); row.names(M.SSASV.abs.nr.t)="M.SSASV.abs.nr"    # transposed to fit combined table below
  M.SSASV.abs.pc.t=t(M.SSASV.abs.pc);row.names(M.SSASV.abs.pc.t)="M.SSASV.abs.pc"     # transposed to fit combined table below
  if (SaveWholeworkspace=='N') {rm(M.SSASV.abs,M.SSASV.abs.nr,M.SSASV.abs.name,M.SSASV.abs.pc)}
} else if (ncol(M.SSASV.abs>1)) {
  M.SSASV.abs.nr=rowSums(M.SSASV.abs)                                                 # Number of SSASVabs per sample
  M.SSASV.abs.pc=round(M.SSASV.abs.nr/M.sobs*100,1)                                   # % SSASVabs
  M.SSASV.abs.nr.t=t(M.SSASV.abs.nr); row.names(M.SSASV.abs.nr.t)="M.SSASV.abs.nr"    # transposed to fit combined table below
  M.SSASV.abs.pc.t=t(M.SSASV.abs.pc);row.names(M.SSASV.abs.pc.t)="M.SSASV.abs.pc"     # transposed to fit combined table below
  if (SaveWholeworkspace=='N') {rm(M.SSASV.abs,M.SSASV.abs.nr,M.SSASV.abs.pc)}
}

# #--- 4.1.6. Number of Relative Singletons ----------------------------------------------------------------------------------------------------------------------------------------------

M.no.SSASVs.abs=M[,-(which(apply(M,2,sum)==1))]                                                       # removes absolute singletons
M.SSASVs.rel=as.matrix(M.no.SSASVs.abs[,which(apply(M.no.SSASVs.abs,2,function(x) any(x==1))==TRUE)]) # calculates relative Single Sequence ASVs (SSASVrel), i.e. relative singletons. (ASVs that occur once in one sample and at least once in at least one other sample)

if (is.null(ncol(M.SSASVs.rel))){                                                       # if dataset has no SSASVrel
  M.SSASVs.rel.nr=rep(0,NR); names(M.SSASVs.rel.nr)=names(M.sobs)
  M.SSASVs.rel.pc=rep(0,NR); names(M.SSASVs.rel.pc)=names(M.sobs)
  M.SSASVs.rel.nr.t=t(M.SSASVs.rel.nr); row.names(M.SSASVs.rel.nr.t)="M.SSASVs.rel.nr"  # transposed to fit combined table below
  M.SSASVs.rel.pc.t=t(M.SSASVs.rel.pc); row.names(M.SSASVs.rel.pc.t)="M.SSASVs.rel.pc"  # transposed to fit combined table below
  if (SaveWholeworkspace=='N') {rm(M.SSASVs.rel.nr,M.SSASVs.rel.pc)}
} else if (ncol(M.SSASVs.rel)==1) {                                                     # if dataset has exactly 1 SSASVrel
  M.SSASVs.rel.name=apply(M.no.SSASVs.abs, 2, function(x) any(x==1)==TRUE)              # Gives the names of all ASVs. TRUE or FALSE indicate if they are singeltons or not.
  M.SSASVs.rel.name=which(M.SSASVs.rel.name, arr.ind = FALSE, useNames = TRUE)          # only keeps SSASVrel
  colnames(M.SSASVs.rel)=as.character(names(M.SSASVs.rel.name))                         # renames the ASVs
  M.SSASVs.rel.nr = matrix(NA, nrow=NR, ncol=1)                                         # matrix to be filled by the function, NR:number of samples
  for (j in 1:NR) {                                                                     # counts the number of SSASVrel per sample and fills the matrix
    M.SSASVs.rel.nr[j,]= sum(M.SSASVs.rel[j,]==1)
  }  
  M.SSASVs.rel.pc=round(M.SSASVs.rel.nr/M.sobs*100,1)                                   # % SSASVrel
  M.SSASVs.rel.nr.t=t(M.SSASVs.rel.nr); row.names(M.SSASVs.rel.nr.t)="M.SSASVs.rel.nr"  
  M.SSASVs.rel.pc.t=t(M.SSASVs.rel.pc); row.names(M.SSASVs.rel.pc.t)="M.SSASVs.rel.pc"  
  if (SaveWholeworkspace=='N') {rm(M.SSASVs.rel.nr,M.SSASVs.rel.pc,M.SSASVs.rel.name)}
} else if (ncol(M.SSASVs.rel)>1) {                                                      # if dataset has more than 1 SSASVrel
  M.SSASVs.rel.nr = matrix(NA, nrow=NR, ncol=1)                                         
  for (j in 1:NR) {                                                                    
    M.SSASVs.rel.nr[j,]= sum(M.SSASVs.rel[j,]==1)                                       
  }    
  M.SSASVs.rel.pc=round(M.SSASVs.rel.nr/M.sobs*100,1)                                   
  M.SSASVs.rel.nr.t=t(M.SSASVs.rel.nr); row.names(M.SSASVs.rel.nr.t)="M.SSASVs.rel.nr"  
  M.SSASVs.rel.pc.t=t(M.SSASVs.rel.pc); row.names(M.SSASVs.rel.pc.t)="M.SSASVs.rel.pc"  
  if (SaveWholeworkspace=='N') {rm(M.SSASVs.rel.nr,M.SSASVs.rel.pc)}
}

# #--- 4.1.7. Combines diversity indices in one table ---------------------------------------------------------------------------------------------------------------------------------

# # creates table of averaged indices, rounds and reformats table
M.diversity=rbind(M.ASVs,M.sobs.r.mean, M.chao1.r.mean, M.invs.r.mean, M.SSASV.abs.nr.t, M.SSASVs.rel.nr.t, M.shan.r.mean,M.SSASV.abs.pc.t, M.SSASVs.rel.pc.t)  # comes from M.seq.tax.subset.ord.ASV; is not ordered based on Grouping1. 
M.diversity.reord=t(rbind(M.diversity[1,],M.sobs,M.diversity[5,],M.diversity[8,],M.diversity[6,],M.diversity[9,],rep(M.minreads,ncol(M.diversity)),M.diversity[2:4,],M.diversity[7,]))
M.diversity.reord[,c(1:3,5,7:9)]=round(M.diversity.reord[,c(1:3,5,7:9)],digits=0)                                                                                   # rounds indices
M.diversity.reord[,c(4,6,10,11)]=round(M.diversity.reord[,c(4,6,10,11)],digits=1)                                                                                       # rounds indices
colnames(M.diversity.reord)=c("Total Reads (unsubsampled)",'Observed ASVs - Richness (unsubsampled)',"Absolute Single Sequence ASVs (unsubsampled)","Absolute Single Sequence ASVs in % (unsubsampled)","Relative Single Sequence ASVs (unsubsampled)","Relative Single Sequence ASVs in % (unsubsampled)",'Total Reads (subsampled)',"Observed ASVs - Richness (subsampled)","Estimated ASVs - Chao1 (subsampled)","Inverse Simpson Diversity (subsampled)","Shannon Entropy (subsampled)")
write.table(M.diversity.reord,file.path(PathToVisuaRAnalysis,"Alpha_Diversity",paste(VisuaRProjectName,"_DiversityIndices.txt",sep="")),sep='\t',col.names = NA)

if (SaveWholeworkspace=='N') {rm(M.ASVs,M.SSASV.abs.nr.t,M.SSASV.abs.pc.t,M.SSASVs.rel.nr.t,M.SSASVs.rel.pc.t,M.diversity.reord)}

# #===  4.2. Visualizes diversity indices ================================================================================================================================================

cat('\n\n4.2. Visualizes diversity indices',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# #===  4.2.1. Boxplots of diversity indices ================================================================================================================================================

cat('\n\n4.2.1. Boxplots of diversity indices',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)


t.M.diversity=t(M.diversity)                                                                                    # Transposes diversity table
if (SaveWholeworkspace=='N') {rm(M.diversity)}
t.M.diversity.groups=cbind(M.metadata.subset.noNA,t.M.diversity)                                                # adds group vectors to diversity index table. 
if (SaveWholeworkspace=='N') {rm(t.M.diversity)}
t.M.diversity.ord=t.M.diversity.groups[order(t.M.diversity.groups[,MapCol1], t.M.diversity.groups[,MapCol2]),]  # orders table using Grouping1 and Grouping2
if (SaveWholeworkspace=='N') {rm(t.M.diversity.groups)}

SOBS.plot=ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1], y=t.M.diversity.ord$M.sobs.r.mean))      # plots Observed ASVs, using MapCol1 (Grouping1) as grouping
SOBS.plot=SOBS.plot+
  stat_boxplot(geom="errorbar",width=0.1)+ 
  geom_boxplot(fill=M.col,outlier.size=1) +                                                                     # M.col is a list of as many colours as present in Grouping1. The first colour is assigned to the first Grouping in alphabetical/increasing order. 
  labs(x=NULL,y="Observed ASVs") + 
  stat_summary(fun.y=mean,colour="black",shape=17,geom="point") +                                               # adds a mean triangle to each boxplot
  stat_n_text()+                                                                                                # adds the number of samples per group at the bottom of the plot
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5),                             
        panel.background = element_blank())
if (length(M.projects.unique)>1){                                                                               # adds significance bars above the boxplots for each combination of Grouping1
  SOBS.plot=SOBS.plot+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.6,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
print(SOBS.plot)

Shan.plot=ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1], y=(t.M.diversity.ord$M.shan.r.mean)))    # plots Shannon
Shan.plot=Shan.plot+
  stat_boxplot(geom="errorbar",width=0.1)+ 
  geom_boxplot(fill=M.col,outlier.size=1) +
  labs(x=NULL, y="Shannon Entropy") + 
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17)+
  stat_n_text()+
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5),
        panel.background = element_blank())
if (length(M.projects.unique)>1){ 
  Shan.plot=Shan.plot+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.6,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
print(Shan.plot)

Simp.plot=ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1], y=(t.M.diversity.ord$M.invs.r.mean)))    # plots Inverse Simpson
Simp.plot=Simp.plot+
  stat_boxplot(geom="errorbar",width=0.1)+ 
  geom_boxplot(fill=M.col,outlier.size=1) +
  labs(x=NULL, y="Inverse Simpson Diversity", col=Grouping1) + 
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17) +
  stat_n_text()+
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())
if (length(M.projects.unique)>1){
  Simp.plot=Simp.plot+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.6,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
print(Simp.plot)

Chao1.plot=ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1], y=(t.M.diversity.ord$M.chao1.r.mean)))    # plots Chao1 richness
Chao1.plot=Chao1.plot+
  stat_boxplot(geom="errorbar",width=0.1)+ 
  geom_boxplot(fill=M.col,outlier.size=1) +
  labs(x=NULL, y="Chao1 Richness", col=Grouping1) + 
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17) +
  stat_n_text()+
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())
if (length(M.projects.unique)>1){
  Chao1.plot=Chao1.plot+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.6,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
print(Chao1.plot)

pdf(file.path(CompositionPlotsPath,paste(VisuaRProjectName,"_DiversityIndices.pdf")),height=4,width=length(M.projects.unique)*3,useDingbats=F)
Diversity.Multiplot=ggarrange(SOBS.plot,Shan.plot,Simp.plot,Chao1.plot,ncol=4,nrow=1,align='hv')                                                 # creates a multiplot of the diversity indices
annotate_figure(Diversity.Multiplot,top = text_grob(paste(VisuaRProjectName," Alpha Diversity"), color = "black", size = 14))
print(Diversity.Multiplot)
dev.off()

if (SaveWholeworkspace=='N') {rm(SOBS.plot,Shan.plot,Simp.plot,Chao1.plot,Diversity.Multiplot)}
closeAllConnections() # closes all currently open connections.

# #===  4.2.2. Violinplots of diversity indices ================================================================================================================================================

cat('\n\n4.2.2. Violinplots of diversity indices',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

SOBS.plot.vio = ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1],y=t.M.diversity.ord$M.sobs.r.mean,fill=t.M.diversity.ord[,MapCol1])) 
if (length(M.projects.unique)>1){ 
  SOBS.plot.vio=SOBS.plot.vio+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.1,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
SOBS.plot.vio=SOBS.plot.vio + 
  geom_violin(trim=T) + 
  scale_fill_manual(values=M.col)+ 
  geom_dotplot(binaxis='y',
               stackdir='center',
               dotsize =0.5,
               fill="black",
               method='histodot',
               binwidth = ((as.numeric(range(t.M.diversity.ord$M.sobs.r.mean)[2]-range(t.M.diversity.ord$M.sobs.r.mean)[1]))/40),
  ) +
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17)+
  theme(legend.position = "none")+
  labs(x=NULL,y="Observed ASVs") + 
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())

print(SOBS.plot.vio)

Shan.plot.vio = ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1],y=t.M.diversity.ord$M.shan.r.mean,fill=t.M.diversity.ord[,MapCol1])) 
if (length(M.projects.unique)>1){ 
  Shan.plot.vio=Shan.plot.vio+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.1,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
Shan.plot.vio=Shan.plot.vio + 
  geom_violin(trim=T) + 
  scale_fill_manual(values=M.col)+ 
  geom_dotplot(binaxis='y',
               stackdir='center',
               dotsize =0.5,
               fill="black",
               method='histodot',
               binwidth = ((as.numeric(range(t.M.diversity.ord$M.shan.r.mean)[2]-range(t.M.diversity.ord$M.shan.r.mean)[1]))/40),
  ) +
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17)+
  theme(legend.position = "none")+
  labs(x=NULL,y="Shannon Entropy") + 
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())
print(Shan.plot.vio)

Simp.plot.vio = ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1],y=t.M.diversity.ord$M.invs.r.mean,fill=t.M.diversity.ord[,MapCol1])) 
if (length(M.projects.unique)>1){ 
  Simp.plot.vio=Simp.plot.vio+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.1,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
Simp.plot.vio=Simp.plot.vio + 
  geom_violin(trim=T) + 
  scale_fill_manual(values=M.col)+ #fills the boxplot
  geom_dotplot(binaxis='y',
               stackdir='center',
               dotsize =0.5,
               fill="black",
               method='histodot',
               binwidth = ((as.numeric(range(t.M.diversity.ord$M.invs.r.mean)[2]-range(t.M.diversity.ord$M.invs.r.mean)[1]))/40),
  ) +
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17)+
  theme(legend.position = "none")+
  labs(x=NULL,y="Inverse Simpson Diversity") + 
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())
print(Simp.plot.vio)



Chao1.plot.vio = ggplot(t.M.diversity.ord, aes(x=t.M.diversity.ord[,MapCol1],y=t.M.diversity.ord$M.chao1.r.mean,fill=t.M.diversity.ord[,MapCol1])) 
if (length(M.projects.unique)>1){ 
  Chao1.plot.vio=Chao1.plot.vio+ 
    geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord[,MapCol1]))),2,simplify=F),
                vjust = 0.1,
                step_increase=0.06,
                size=0.2,
                textsize=3,
                tip_length=0.01,
                map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2))}
Chao1.plot.vio=Chao1.plot.vio + 
  geom_violin(trim=T) + 
  scale_fill_manual(values=M.col)+ #fills the boxplot
  geom_dotplot(binaxis='y',
               stackdir='center',
               dotsize =0.5,
               fill="black",
               method='histodot',
               binwidth = ((as.numeric(range(t.M.diversity.ord$M.chao1.r.mean)[2]-range(t.M.diversity.ord$M.chao1.r.mean)[1]))/40),
  ) +
  stat_summary(fun.y=mean,colour="black",geom="point",shape=17)+
  theme(legend.position = "none")+
  labs(x=NULL,y="Chao1 Richness") + 
  theme(axis.title = element_text(colour="black",size=16),
        axis.text = element_text(colour="black"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        axis.line=element_line(colour="black"),
        axis.ticks=element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour="grey",size=0.25),
        panel.grid.minor.y = element_line(),
        panel.border = element_rect(fill=NA,colour = "black",size=0.5), 
        panel.background = element_blank())
print(Chao1.plot.vio)

pdf(file.path(CompositionPlotsPath,paste(VisuaRProjectName,"_DiversityIndices_Violin.pdf")),height=4,width=length(M.projects.unique)*3,useDingbats=F)
Diversity.Multiplot.vio=ggarrange(SOBS.plot.vio,Shan.plot.vio,Simp.plot.vio,Chao1.plot.vio,ncol=4,nrow=1,align='hv')                                  # creates a multiplot of the diversity indices
annotate_figure(Diversity.Multiplot.vio,top = text_grob(paste(VisuaRProjectName," Alpha Diversity"), color = "black", size = 14))
print(Diversity.Multiplot.vio)
dev.off()

if (SaveWholeworkspace=='N') {rm(SOBS.plot.vio,Shan.plot.vio,Simp.plot.vio,Chao1.plot.vio,Diversity.Multiplot.vio)}
if (SaveWholeworkspace=='N') {rm(t.M.diversity.ord)}

closeAllConnections() # closes all currently open connections.

# #======= 4.3. Species accumulation curves ====================================================================================================================================================

cat('\n\n4.3. Species accumulation curves',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

# #------- 4.3.1. Calculates Species accumulation curve using samples -----------------------------------------------------------------------------------------------------------------

cat('\n4.3.1. Calculates Species accumulation curve using samples',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.spec=specaccum(M, method="random", permutations=NI, conditioned=T, gamma="Chao")
pdf(file.path(CompositionPlotsPath,paste(VisuaRProjectName,"_Species_Accumulation.pdf")),useDingbats=F)
par(mar=c(5,2,5,7), xpd=TRUE)                                                                             # mar: margin sizes in the following order: bottom, left, top, and right.
plot(M.spec,ci.type="poly", col="grey", lwd=3, ci.lty=0, ci.col=alpha('grey',alpha=0.5),main=c(VisuaRProjectName,'Species Accumulation_Total'),ylab='Species richness',xlab='Samples')
boxplot(M.spec,col="grey", add=TRUE, outpch=19,outcex=0.35,boxwex=0.25,medlwd=1.8,whisklty=1)             # medlwd: width of the median line, boxwex:a scale factor applied to all boxes (to make them narrower or wider), whisklty=1 sets whisker line type to line (not dashed as default), outcex: size of outlier points, outpch=19: outlier points are filled dots
par(new=T)
plot.new()
legend(x=1.05,y=1.0,legend=c('All samples'),fill=c('grey'),bty ='n',cex=0.75,text.font = 2)
dev.off()

# #------- 4.3.2. Calculate Species accumulation curve using Grouping1 -----------------------------------------------------------------------------------------------------------------------

cat('\n4.3.2. Calculate Species accumulation curve using ',Grouping1,' (Grouping1).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.plus.metadata=cbind(as.data.frame(M),M.metadata.subset.noNA[,Grouping1]) # adds Grouping1 column to M

pdf(file.path(CompositionPlotsPath,paste(VisuaRProjectName,"_Species_Accumulation_Groups.pdf")),useDingbats=F)
par(mar=c(5,2,5,7), xpd=TRUE) 
plot(M.spec,ci.type="poly", col='grey', lwd=3, ci.lty=0, ci.col=alpha('grey',alpha=0.2),main=c(VisuaRProjectName,'Species Accumulation_Total_and_Grouped'),ylab='Species richness',xlab='Samples')  # The accumulation curve by sample is plotted first
for (i in 1:length(M.projects.unique)) {                                                                                                                                                            # adds accumulation curves by Group as in Grouping1 to the above curve
  SubsetDataframe=M.plus.metadata[grep(M.projects.unique.ord[i],M.plus.metadata[,ncol(M.plus.metadata)],invert=F),]                                                                                 # only takes samples belonging to the first grouping (alphabetically)
  SubsetDataframe=SubsetDataframe[,-ncol(SubsetDataframe)]                                                                                                                                          # deletes the column with the Grouping1 information
  loopvariable=paste('Specaccum.Group.',M.projects.unique.ord[i],sep='')                                                                                                                            # creates a new 'loopvariable' for each round
  to.specaccum.plot=assign(loopvariable,specaccum(SubsetDataframe, method="random", permutations=NI, conditioned=T, gamma="Chao"))                                                                  # assign: assigns the value calculated in specaccum (the respective species accumulation of the group calculating at this moment) to the variable saved in loopvariable (not to loopvariable itself). In addition the result is saved in the variable 'to.specaccum.plot'. 
  plot(to.specaccum.plot,ci.type="poly", col=M.col[i], lwd=3, ci.lty=0, ci.col=alpha(M.col[i],alpha=0.25),add=T)                                                                                    # adds the specaccum plot of each group to the plot
  rm(to.specaccum.plot)
}
par(new=T)
plot.new()
legend(x=1.05,y=1.0,legend=c('All',M.projects.unique.ord),fill=c('grey',M.col),bty ='n',cex=0.75,text.font = 2) # adds a legend to the above created plot
dev.off()

if (SaveWholeworkspace=='N') {rm(M.spec,M.plus.metadata,SubsetDataframe,loopvariable)}

closeAllConnections() # closes all currently open connections.

# #================== 5. Beta Diversity =================================================================================================================================================

cat('\n\n5. Beta Diversity',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.grouped=rowsum(M, group=M.groups, reorder = TRUE)   # sums up ASV abundances of a condition using Grouping1 (columns:ASVs, rows:groups,ordered alphabetically/increasing, fill: ASV reads)
M.grouped.pa=decostand(M.grouped, method='pa')        # creates presence absence table of the Group  by ASV table


# #========== 5.1. Anaylsis of Similarities (ANOSIM) and calculation of Non-metric Multi-Dimensional Scaling (NMDS) ordinations. =====================================================

cat('\n5.1.  Anaylsis of Similarities (ANOSIM) and calculation of Non-metric Multi-Dimensional Scaling (NMDS) ordinations.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

NMDSPlotsPath=file.path(PathToVisuaRAnalysis,"Beta_Diversity","NMDSPlots")
dir.create(file.path(PathToVisuaRAnalysis,"Beta_Diversity","NMDSPlots"))

# #------- 5.1.1. ANOSIM ----------------------------------------------------------------------------------------------------------------------------------------------------------------

if(length(M.projects.unique)>1){
  M.ano=anosim(M.rel, permutations=999, grouping=M.groups, distance="bray")                                                           # Anosim calculation of different groups . M.rel is a sample by ASV table filled with the relative abundances. The samples are not in the order of the sorted Grouping1. M.groups is in the same 'unordered' format
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Anosim.pdf")),useDingbats=F)
  plot(M.ano,col=c("black",M.col),xlab=NULL,ylab='Dissimilarity Rank',main=paste('Anosim\n',VisuaRProjectName),cex.main=0.8,xaxt='n') # plots the result of ANOSIM
  axis(1, at=1:(length(M.projects.unique)+1),labels=c('Between',M.projects.unique.ord))                                               # renames the axis with the Grouping1 names
  dev.off()
  cat('\n6.5.1.1. Analysis of Similarities (ANOSIM)',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat("\nAnosim_Dissimilarity: ",M.ano$dissimilarity,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat("\nAnosim_statistic_R: ",M.ano$statistic,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat("\nAnosim_Significance: ",M.ano$signif,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
  cat("\nAnosim_Permutations: ",M.ano$permutations,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
}

if (SaveWholeworkspace=='N') {rm(M.rel)}

# #--- 5.1.2. NMDS Calculation And Visualization ----------------------------------------------------------------------------------------------------------------------------------------
cat('\n\n5.1.2. NMDS calculation and Visualization',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nNote: Using metaMDS with a distance matrix is the same as using monoMDS without a distance matrix. You create a distance matrix using vegdist(method=bray) and apply metaMDS to this distrance matrix.' ,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.dist=vegdist(M, method="bray") # Calculates distance matrix. Samples in M are not ordered based on Grouping1

# # Community Similarity in %
Max.Com.Sim=(1-(min(M.dist)))*100 
Min.Com.Sim=(1-(max(M.dist)))*100
Mean.Com.Sim=(1-(mean(M.dist)))*100
cat('\nNote: Maximum, Mean and Minimal Community Similarity in % (as shown in the NMDS plots) refers to the comparison of any two samples using the bray-curtis distance matrix.' ,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)

M.mMDS=metaMDS(M.dist) # takes distances from n-dimensional M.dist (n=number of samples-1) and creates a 2-dimensional ordination. 
if (SaveWholeworkspace=='N') {rm(M.dist)}
saveRDS(M.mMDS,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_distanceMatrix.rds",sep="")))
cat("\nNote: In NMDS ordinations the distance between two samples represents the distance of their underlying communities.\nNMDS distances are relative measures and thus, do not need an axis.\nAxes in NMDS ordinations are only used when NMDS is combined (i.e. superimposed on the NMDS ordination) with another method.",file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)

M.mMDS.stress=M.mMDS$stress 
cat('\nNMDS_Stress_Level: ',M.mMDS.stress,file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep='',append=TRUE)
cat('\nNote: Stress values are a measure for the goodness of the NMDS (i.e. how realistic is the n-dimensional space represented in the 2-dimensional plot).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)
cat('\nNote: R and p values shown in the NMDS plots are calculated using ANOSIM. Please consider limitations and interpretation of the ANOSIM test.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)
cat('\nNote: NMDS featuring text of the samples is meant for reference and to look up specific samples yet is not optimized for display.',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)

# # NMDS featuring sample names 
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_Text.pdf",sep="")),height=20,width=20,useDingbats=F)
par(mar=c(3,3,3,3),xpd=T)
ordiplot(M.mMDS, type="t", display="sites",ylab="",xlab="",xaxt='n',yaxt='n')                                     # type='t': shows the sites as text, display='sites: shows samples, y and x axis removed
dev.off()

# # NMDS featuring Groups
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_dots.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=T)                                                                                     # mar: margin sizes in the following order: bottom, left, top, and right. There is more space at the right hand site for the figure legend, xpd='T': all plotting is clipped to the figure region (not only to the plot region). This allows to place figure legens outside of the plot
ordiplot(M.mMDS, type='n', display="sites",ylab='',xlab='',xaxt='n',yaxt='n')                                     # type='n': sites are invisible and can be filled with points command
points (M.mMDS, col=M.colvec, pch=19)                                                                             # shows sites as dots colored according to groups (pch=19: filled dot)
par(new=T)                                                                                                        # only using par(new=T) it is possible to plot the figure legend independent of the axis scaling to a fixed place inside the pdf.
plot.new()
if(length(M.projects.unique)>1){                                                                                  # only if more than one grouping
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                                            # text.font=2: prints text in bold, bty='n': no box will be drawn around the legend. Alternative: bty='o'
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75)                             # cex: changes the size of the Stress value shown in the plot
} else {
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

# # NMDS featuring Shannon Entropy
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_Shannon.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=TRUE)
ordiplot (M.mMDS, display = 'si', type = 'n',ylab='',xlab='',xaxt='n',yaxt='n')
points (M.mMDS, col=M.colvec, pch=19, cex=M.shan.r.mean/(max(M.shan.r.mean)/2.7))                                 # cex= the point size of each sample represents Shannon Entropy.
par(new=T)                                                                                                        # cex sizes are scaled to fit the plot dimensions using a scaling factor based on the maximum average Shannon Entropy and a constant (2.7)
plot.new()                                                                                                        # circle sizes are meant to visualize relative differences in Shannon. Not absolute values. 
if(length(M.projects.unique)>1){ 
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2) 
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Shannon Entropy',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75)
} else {
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Shannon Entropy',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

if (SaveWholeworkspace=='N') {rm(M.shan.r.mean)}

# # NMDS featuring ASV Richness
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_Richness.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=TRUE)
ordiplot (M.mMDS, display = 'si', type = 'n',ylab='',xlab='',xaxt='n',yaxt='n')
points (M.mMDS, col=M.colvec, pch=19, cex=M.sobs.r.mean/(max(M.sobs.r.mean)/2.7)) 
par(new=T)
plot.new()
if(length(M.projects.unique)>1){ # only if more than one grouping
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='ASV Richness',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75)
} else {
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='ASV Richness',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

if (SaveWholeworkspace=='N') {rm(M.sobs.r.mean)}

# # NMDS featuring Inverse Simpson Diversity
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_InvSimpson.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=TRUE)
ordiplot (M.mMDS, display = 'si', type = 'n',ylab="",xlab="",xaxt='n',yaxt='n')
points (M.mMDS, col=M.colvec, pch=19, cex=M.invs.r.mean/(max(M.invs.r.mean)/3.5)) 
par(new=T)
plot.new()
if(length(M.projects.unique)>1){ 
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2) 
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Inverse Simpson',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
} else {
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Inverse Simpson',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

if (SaveWholeworkspace=='N') {rm(M.invs.r.mean)}

# # NMDS featuring Chao1 Richness 
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_Chao1.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=TRUE) 
ordiplot (M.mMDS, display = 'si', type = 'n',ylab="",xlab="",xaxt='n',yaxt='n')
points (M.mMDS, col=M.colvec, pch=19, cex=M.chao1.r.mean/(max(M.chao1.r.mean)/3.5))
par(new=T)
plot.new()
if(length(M.projects.unique)>1){ 
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2) 
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Chao1 Richness',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
} else {
  legend(1.01,0.58,legend='Dotsize',bty='n',cex = 0.75,text.font = 2)
  legend(1.01,0.535,legend='Chao1 Richness',bty='n',cex = 0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

if (SaveWholeworkspace=='N') {rm(M.chao1.r,M.chao1.r.mean)}

# # NMDS featuring Ordiellipse&Ordispider 
cat('\nNote: Ordispider: each Sample is connected to the weighted average mean of the within group distances (centroid).',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)
cat('\nNote: Ordiellipse: visualizes one standard deviation from the groups centroids. ',file=(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_Analysis_log.txt',sep=''))),sep="\t",append=TRUE)

pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_spider_ellipse.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=F)                                                                                             # first, xpd=F because if not the ellipses might be drawn outside of the plot region
ordiplot (M.mMDS, display = 'si', type = 'n',ylab="",xlab="",xaxt='n',yaxt='n')
for (i in seq (1, NR)) ordiellipse (M.mMDS, groups = M.groups, show.groups = i, col = M.col[i], label = F,lwd=1.5,draw='polygon',border='NA',alpha=60,kind='sd')
for (i in seq (1, NR)) ordispider (M.mMDS, groups = M.groups, show.groups = i, col = M.col[i], label = F,lwd=1.5)
points (M.mMDS, col=M.colvec, pch=19,cex=1) 
par(new=T,xpd=T)
plot.new()
if(length(M.projects.unique)>1){
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2) 
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
} else {
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()


# # NMDS featuring Ordihull 
pdf(file.path(NMDSPlotsPath,paste(VisuaRProjectName,"_NMDS_Ordihull.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=TRUE) 
ordiplot (M.mMDS, display = 'si', type = 'n',ylab="",xlab="",xaxt='n',yaxt='n')
points (M.mMDS, col=M.colvec, pch=19,cex=0.5) 
for (i in seq (1, NR)) ordihull (M.mMDS, groups = M.groups, show.groups = i, col = M.col[i], label = F)
par(new=T)
plot.new()
if(length(M.projects.unique)>1){ 
  legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)
  legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.21,legend=paste0('ANOSIM'),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.17,legend=paste0('R: ',round(M.ano$statistic,3),'\np: ',M.ano$signif),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
} else {
  legend(1.01,0.46,legend=paste0("Community\nSimilarity [%]"),bty='n',cex=0.75,text.font = 2)
  legend(1.01,0.4,legend=paste0('max: ',round(Max.Com.Sim,1),'\nmean: ',round(Mean.Com.Sim,1),'\nmin: ',round(Min.Com.Sim,1)),bty='n',cex=0.75)
  legend(1.01,0.03,legend=paste0("Stress: ",round(M.mMDS.stress,3)),bty='n',cex=0.75) 
}
dev.off()

if(length(M.projects.unique)>1){
  if (SaveWholeworkspace=='N') {rm(M.ano)}
}

if (SaveWholeworkspace=='N') {rm(M.mMDS)}

closeAllConnections() # closes all currently open connections.

# #------- 5.2. Calculate and visualize grouped rarefaction -----------------------------------------------------------------------------------------------------------------------------

MRR=M.grouped #sample-wise: M; group-wise: M.grouped

M.rar <- rrarefy(MRR,sample=M.minreads) # data frame with estimated species richness using random subsamples of size M.minreads

pdf(file.path(CompositionPlotsPath,paste(VisuaRProjectName,"_RarefactionCurve.pdf",sep="")),height=5,width=6,useDingbats=F)
par(mar=c(5, 5, 5, 7), xpd=F) 
M.rac <- rarecurve(MRR,sample=M.minreads,step=1000,se=T,col=M.col,lwd=3,main=c(VisuaRProjectName,'Rarefaction_Grouped'),cex.main=0.8,label=F) # rarefaction curve for each group. yields richness of all samples at M.minreads
par(new=T,xpd=T)
plot.new()
legend(x=1.05,y=1.05,legend=c(M.projects.unique.ord),fill=c(M.col),bty ='n',cex=0.75,text.font = 2) 
dev.off()

# M.ras <- rareslope(MRR,sample=M.minreads-1) # calculates the slope of rarecurve for each group

if (SaveWholeworkspace=='N') {rm(MRR,M.rar)}

# #------- 5.3. Calculates and visualizes Cluster Dendrograms by Groups ------------------------------------------------------------------------------------------------------------------
# # Default metric is bray curtis distance and average linkage hierarchical clustering.
dir.create(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms'))

if(length(M.projects.unique)>2){  
  M.grouped.dist=vegdist(M.grouped, method="bray")                              # Calculates distance matrix
  M.grouped.clust <- hclust(M.grouped.dist, method="average")                   # hierarchical cluster analysis on a set of dissimilarities
  M.grouped.clust$labels=c(M.projects.unique.ord)                               # renames the groups (names are removed while clustering)
  
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_ClusterDendrogram_grouped.pdf",sep="")),height=5,width=5,useDingbats=F)
  M.grouped.cd=plot(M.grouped.clust,main=paste('Cluster Dendrogram\n',VisuaRProjectName),cex.main=0.8,xlab =NA,ylab = 'Distance' ,sub=NA)
  dev.off()
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_ClusterDendrogram2_grouped.pdf",sep="")),height=5,width=5,useDingbats=F)
  M.grouped.phylo=plot(as.phylo(M.grouped.clust),label.offset = 0.01,cex=0.6,tip.color = c(M.col))
  dev.off()
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_Cluster_fan_grouped.pdf",sep="")),height=5,width=5,useDingbats=F)
  M.grouped.phylo=plot(as.phylo(M.grouped.clust),type='fan',cex=0.4,label.offset = 0.02,tip.color = c(M.col))
  dev.off()
  if (SaveWholeworkspace=='N') {rm(M.grouped.dist,M.grouped.clust,M.grouped.cd,M.grouped.phylo)}
}

if (SaveWholeworkspace=='N') {rm(M.grouped)}

# #------- 5.4. Calculates and visualizes Cluster Dendrograms by samples -----------------------------------------------------------------------------------------------------------------
# # Creates different cluster diagrams of all samples colored in the selected Grouping1
M.sample.dist=vegdist(M, method="bray") 
M.sample.clust <- hclust(M.sample.dist, method="average")

pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_ClusterDendrogram_samples.pdf",sep="")),height=5,width=5,useDingbats=F)
M.sample.cd=plot(M.sample.clust,main=paste('Cluster Dendrogram\n',VisuaRProjectName),cex.main=0.8,xlab =NA,ylab = 'Distance' ,sub=NA,hang=-1)
dev.off()
pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_ClusterDendrogram3_samples.pdf",sep="")),height=5,width=5,useDingbats=F)
M.sample.phylo=plot(as.phylo(M.sample.clust),no.margin = T,cex=0.6,tip.color = c(M.colvec))
dev.off()
pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_Cluster_fan_samples.pdf",sep="")),height=5,width=5,useDingbats=F)
M.sample.phylo=plot(as.phylo(M.sample.clust),type='fan',no.margin = T,cex=0.4,tip.color = c(M.colvec))
dev.off()
pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity','Dendrograms',paste(VisuaRProjectName,"_Cluster_radial_samples.pdf",sep="")),height=5,width=5,useDingbats=F)
M.sample.phylo=plot(as.phylo(M.sample.clust),type='radial',no.margin = T,cex=0.6,tip.color = c(M.colvec))
dev.off()

if (SaveWholeworkspace=='N') {rm(M.sample.dist,M.sample.clust,M.sample.phylo)}

# #------- 5.5. Create UpsetR turnover plot ---------------------------------------------------------------------------------------------------------------------------------------------
# # UpsetR is a tool to visualize overlaps in datasets, it is analogous to Venn diagrams, yet can plot more dimensions
# # Check Conway et al. 2017, doi: 10.1093/bioinformatics/btx364
# # It can be used with samples (M.pa) or groups (M.grouped.pa). For samplewise analysis nsets has to be changed to the number of samples
# # If only specific samples should be visualize use sets=c('name1','name2',...) and set nsets to the number of samples
# # nintersects: number of intersects to show, set to NA to show all 

if(length(M.projects.unique)>1){
  M.up=as.data.frame(t(M.grouped.pa)) 
  colnames(M.up)=M.projects.unique.ord
  M.upset=upset(M.up,mb.ratio = c(0.7, 0.3),nsets = length(M.projects.unique),sets=rev.default(M.projects.unique.ord),keep.order = T, nintersects = NA, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),sets.bar.color = rev.default(M.col),empty.intersections = "on",point.size = 5,text.scale=2) # keep.order=T: Keep Groups in the order entered using the sets parameter. The default is FALSE, which orders the sets by their sizes.,mb.ratio: Ratio between matrix plot and main bar plot
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_UpsetPlot.pdf",sep="")),height=6,width=15,useDingbats=F)
  M.upset
  print(M.upset)
  dev.off()
  if (SaveWholeworkspace=='N') {rm(M.upset)}
}

# #------- 5.6. Create Venn Diagram -----------------------------------------------------------------------------------------------------------------------------------------------------

if(length(M.projects.unique)>1){
  pdf(file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_VennDiagram.pdf",sep="")),height=5,width=5,useDingbats=F)
  venn(M.up,zcolor = M.col,opacity = 0.4,ilcs=0.85,sncs=1)
  dev.off()
  if (SaveWholeworkspace=='N') {rm(M.up)}
}

# #------- 5.7. Pairwise comparison of groups -------------------------------------------------------------------------------------------------------------------------------------------

if(length(M.projects.unique)>1){
  M.compare=M.grouped.pa  # group-wise analysis, M.grouped.pa is a presence absence matrix from grouped matrix of M
} else {
  M.compare=M.pa          # Sample-wise analysis
}

if (SaveWholeworkspace=='N') {rm(M.pa,M.grouped.pa)}

# # number of unique ASVs of GroupX when comapred with GroupY... 
unique_for_group=matrix(NA, (nrow(M.compare)^2), 3) # creates a matrix of (number of Groups)^2 rows and 3 columns

# # number of total unique ASVs 
unique=matrix(NA,nrow(M.compare),nrow(M.compare))
rownames(unique)= M.projects.unique.ord
colnames(unique) = M.projects.unique.ord

# # number of shared OTUs 
common=matrix(NA,nrow(M.compare),nrow(M.compare)) 
rownames(common)=M.projects.unique.ord
colnames(common)=M.projects.unique.ord

# # number of total ASVs 
total=matrix(NA,nrow(M.compare),nrow(M.compare)) 
rownames(total)= M.projects.unique.ord
colnames(total) = M.projects.unique.ord

# # Percentage of shared ASVs ((shared number of ASVs/total number of ASVs)* 100)
percentage=matrix(NA,nrow(M.compare),nrow(M.compare)) 
rownames(percentage)= M.projects.unique.ord
colnames(percentage) = M.projects.unique.ord

b=vector() # creates an empty vector
p=vector()

rowid = 0                                                               # sets rowid to 0. This will be used to identify the row to fill the values into (first round first row, second round second row....)
for (r in (1 : (nrow(M.compare)))) {                                    # number of repeats = number of groups
  for (s in (1 : (nrow(M.compare)))) {                                  # number of repeats = number of groups
    for (c in (1 : (ncol(M.compare)))) {                                # number of repeats = number of ASVs
      b[c]=sum (M.compare[r,c],M.compare[s,c])                          # vector b is filled with the sum of cell [1,1] and cell [1,1] (=1 if only in one, =2 if in both, =0 if in none)
      p[c] = if (M.compare[r,c]==1 & M.compare[s,c]== 0) 1 else 0       # vector p is filled with 1 if the ASV occur in cell 1 but not in cell 2 and with 0 if 0 in both or 1 in both
    }
    rowid = rowid + 1;
    unique_for_group[rowid,1] = paste(M.projects.unique.ord[r])         # prints the name of the current sample[r] in the current row[rowid]
    unique_for_group[rowid,2] = paste(M.projects.unique.ord[s])         # prints the name of the current sample[s] in the current row[rowid]
    unique_for_group[rowid,3] = sum(p)                                  # prints the sum of the unique ASVs for sample [r] (do occur in samle r but not in sample [s])
    
    unique[r,s]=length (which (b==1) )                                  # fills unique table with the number of how many ASVs only occur in sample [r] or sample [s]
    common[r,s]= length (which (b==2) )                                 # fills common table with the number of how many ASVS occur in both samples ([r] and [s])
    total[r,s]= sum (unique[r,s], common[r,s])                          # fills total table with the sum of the unique and common ASVs for each sample
    percentage[r,s]=round(common[r,s]/total[r,s]*100,0)                 # fills percentage table with the percentage of shared ASVs for each sample
  }
}

if (SaveWholeworkspace=='N') {rm(M.compare)}

write.table(total,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Total_ASVs_Groups.txt",sep="")),sep='\t',col.names = NA)
write.table(unique,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Unique_ASVs_Groups.txt",sep="")),sep='\t',col.names = NA)
write.table(common, file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Shared_ASVs_Groups.txt",sep="")),sep='\t',col.names = NA)
write.table(percentage,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Percent_Shared_ASVs_Groups.txt",sep="")),sep='\t',col.names = NA)
write.table(unique_for_group,file.path(PathToVisuaRAnalysis,'Beta_Diversity',paste(VisuaRProjectName,"_Unique_ASVs_per_Grouping.txt",sep="")),sep='\t',col.names = NA) # results of "unique_for_group" should read as: sampleI (first column), when compared to sampleJ (second column), has "X Number" (third column) of unique ASVs

if (SaveWholeworkspace=='N') {rm(total,unique,common,percentage,unique_for_group)}

# #=== 6. Epilogue ============================================================================================================================================================================
save.image(file.path(PathToVisuaRAnalysis,paste(VisuaRProjectName,'_VisuaR','.RData',sep='')))
