basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
package.list <- setdiff(package.list,basic.packages)
if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
rm(list=ls()); closeAllConnections(); gc()
# # https://microbiome.github.io/tutorials/Core.html
# # first create phlyoseq object from VisuaR analysis (VisuaR_addon_createPhyloseq.R)

VisuaRName <- "" # Name of your VisuaR analysis.
VisuaRPath <- "" # Path to your VisuaR analysis.

physeq <- readRDS(paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_phyloseq.rds"))
dir.create(file.path(paste0(VisuaRPath,"/",VisuaRName,"/corebiome")))

Groupcolumn1 <- "" # # e.g. biome
GroupCondition1 <- "" # e.g. Air
Groupcolumn2 <- "" # # e.g. air origin
GroupCondition2 <- "" # # e.g. terrestrial

pathtooutput <- file.path(paste0(VisuaRPath,"/",VisuaRName,"/corebiome/",GroupCondition1,"_",GroupCondition2)) # change the last entry to your liking
dir.create(pathtooutput)


physeq <- phyloseq::subset_samples(physeq, physeq@sam_data[[Groupcolumn1]] == GroupCondition1) # subsets to first grouping
physeq <- phyloseq::subset_samples(physeq, physeq@sam_data[[Groupcolumn2]] == GroupCondition2) # subsets to second grouping


# # keep only taxa with positive sums (necessary if subsetted before)
physeq <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
# # Calculate compositional version of the data (relative abundances)
physeq.rel.orig <- microbiome::transform(physeq, "compositional")
rm(physeq);gc()
tax.levels <- c("Phylum","Class","Order","Family","Genus","Species")


for (i in 0:length(tax.levels)) {
  if (i != 0) {
    physeq.rel <- microbiome::aggregate_taxa(physeq.rel.orig,tax.levels[i])
  } else if (i==0) {
    physeq.rel <- physeq.rel.orig
  }
  gc()
  # # full phyloseq object of core microbiota and then extract taxonomy from phyloseq object
  # detection: average relative abundance at least this high 
  # prevalence: we use include.lowest = T, i.e. if prevalence is 0.5, clades which are 50% prevalent will still be included
    # # calculates community core abundance index per sample, i.e. relative abundance of the core community in each sample
  # pseq.core.abund <- data.table::as.data.table(microbiome::core_abundance(physeq.rel, detection = detection, prevalence = prevalence, include.lowest = T))
  pseq.core.abund <- as.data.frame(microbiome::core_abundance(physeq.rel, detection = 0.001, prevalence = 0.5, include.lowest = T))
  # # returns members of core microbiome
  pseq.core <- data.table::as.data.table(microbiome::core_members(physeq.rel, detection = 0.001, prevalence = 0.5, include.lowest = T))
  
  if (i==0) {
    data.table::fwrite(pseq.core.abund,file=paste0(pathtooutput,"/",VisuaRName,"_core-asv-relabund-persample.csv",sep=""),row.names = T)
    data.table::fwrite(pseq.core,file=(paste0(pathtooutput,"/",VisuaRName,"_core-asv-names.csv",sep="")))
    
  } else {
    data.table::fwrite(pseq.core.abund,file=paste0(pathtooutput,"/",VisuaRName,"_core-",tax.levels[i],"-relabund-persample.csv",sep=""),row.names = T)
    data.table::fwrite(pseq.core,file=(paste0(pathtooutput,"/",VisuaRName,"_core-",tax.levels[i],"-names.csv",sep="")))
  }
  
  rm(pseq.core.abund); gc()
  # Core line plots
  # Determine core microbiota across various abundance/prevalence thresholds with the blanket analysis (Salonen et al. CMI, 2012) based on various signal and prevalences.
  # With compositional (relative) abundances
  det <- c(0,0.001,0.01,0.1, 0.5, 2, 5, 20, 40)/100
  prevalences <- seq(.05, 1, .05)
  core.microbiota.explore <- microbiome::plot_core(physeq.rel, prevalences = prevalences,detections = det, plot.type = "lineplot") +
    ggplot2::xlab("Relative Abundance (%)") +
    ggplot2::theme_bw() 
  if(i==0) {
    core.microbiota.explore <- core.microbiota.explore +
      ggplot2::ylab("Core size (N) - ASVs")
  } else {
    core.microbiota.explore <- core.microbiota.explore +
      ggplot2::ylab(paste0("Core size (N) - ",tax.levels[i],sep=""))
    
  }
  core.microbiota.explore
  if (i==0) {
    pdf(file.path(paste0(pathtooutput,"/",VisuaRName,"_ASVs-coresize-by-relabund.pdf",sep="")),height=5,width=8,useDingbats=F); print(core.microbiota.explore); dev.off()
      } else {
    pdf(file.path(paste0(pathtooutput,"/",VisuaRName,"_",tax.levels[i],"-coresize-by-relabund.pdf",sep="")),height=5,width=8,useDingbats=F); print(core.microbiota.explore); dev.off()
    
  }
  rm(core.microbiota.explore,prevalences,det);gc()
  
  # # Heat map with Prevalence as color code and detection threshold (relative abundance at x axis)
  prevalences <- seq(.05, 1, .05)
  # detections <- round(10^seq(log10(1e-5), log10(.4), length = 10), 3) # 0.05 to 
  detections <- c(0.001,0.01,0.1, 0.5, 2, 5, 20, 40)/100
  
  p1 <- microbiome::plot_core(physeq.rel,
                              plot.type = "heatmap",
                              # colours = gray,
                              prevalences = prevalences,
                              detections = detections, min.prevalence = .2) +
    ggplot2::xlab("Detection Threshold (Relative Sequence Abundance (%))") +
    ggplot2::theme_bw() + 
    ggplot2::ylab(tax.levels[i]) +
    ggplot2::scale_fill_gradientn(colors=colorspace::sequential_hcl(n=7,palette="Plasma"))
  p1
  
  if (i==0) {
    pdf(file.path(paste0(pathtooutput,"/",VisuaRName,"_ASVs-coreheatmap.pdf",sep="")),height=(1+(nrow(pseq.core)*0.4)),width=4,useDingbats=F); print(p1); dev.off()
    
  } else {
    pdf(file.path(paste0(pathtooutput,"/",VisuaRName,"_",tax.levels[i],"-coreheatmap.pdf",sep="")),height=(1+(nrow(pseq.core)*0.4)),width=(4+i),useDingbats=F); print(p1); dev.off()
  }
  
  rm(p1,physeq.rel,detections,prevalences);gc()
}
