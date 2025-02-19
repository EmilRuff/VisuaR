basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
package.list <- setdiff(package.list,basic.packages)
if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
rm(list=ls()); closeAllConnections(); gc()

library(dplyr)
# Install and load the ggtern and RColorBrewer packages if not already installed
# if (!requireNamespace("ggtern", quietly = TRUE)) {
#   install.packages("ggtern")
# }
# if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
#   install.packages("RColorBrewer")
# }

VisuaRName <- "" # enter the name of the VisuaR analysis, e.g., "ProjectName"
VisuaRPath <- "" # enter the path to the VisuaR analysis, including the /VisuaR

numberofcoloredclades <- 5  # number of phyla/genera etc. that should be colored in the ternary plot, there will be the top numbers shown and all others in grey

# # select  the combinations you want to plot. If left like this, it will result in all combinations. Only works if the taxonomy names match the ones below.
# # e.g., levelsOfGrouping "phylum" and lovelsOfDots "Order" creates a ternary plot with one dot per order, colored by the phylum it belongs to.
levelsOfGrouping <- c("kingdom","phylum","class","order","family","genus","species")
levelsOfGroupingCaptial <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species") 
levelsOfDots <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","ASV") 

#########no need to change anything from here on #######################
# # loop through grouping levels
for (i in 1:length(levelsOfGrouping)) {
  levelOfGrouping <- levelsOfGrouping[i]
  levelOfGroupingCaptial <- levelsOfGroupingCaptial[i]
  
  # # loop through dot levels (always starting one level lower than grouping level)
  for (j in i:length(levelsOfDots)) {
    levelOfDots <- levelsOfDots[j]
    
    ASV.by.sample.and.tax <- as.data.frame (data.table::fread(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/",VisuaRName,"_ASVbySample_relabund_final.csv",sep="")))
    ASV.by.sample.and.tax <- ASV.by.sample.and.tax[,-c(1:2)]
    
    # # read in relative abundance of level of interest (this is the level used for coloring/grouping of the ASVs, not in the same sample order as the ASV table
    levelOfGroupingTOP <- as.data.frame(data.table::fread(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/CompositionTables/",VisuaRName,"_",levelOfGrouping,"_relabund.csv",sep="")))
    # # group vector (sorted similar as ASV.by.sample.and.tax, not as levelOfGroupingTOP)
    M.groups <- readRDS(paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_Mgroups.rds",sep=""))
    # # color vector (sorted alphabetically by grouping)
    M.col <- readRDS(paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_Mcol.rds",sep=""))
    # # group names (sorted alphabetically by grouping)
    M.projects.unique.ord <- readRDS(paste0(VisuaRPath,"/",VisuaRName,"/",VisuaRName,"_Mprojectsuniqueord.rds",sep=""))
    
    n.samples <- length(M.groups)
    
    # # create summed table for level of Dots (sums up over e.g. Genera)
    if (levelOfDots != "ASV") {
      # # first we summarise the current level of dots, e.g. Genus together with the grouping column (e.g. phylum), this will reorder the table alphabetically (not the samples!)
      groupingdotassignment <- cbind(ASV.by.sample.and.tax[, levelOfGroupingCaptial],ASV.by.sample.and.tax[, levelOfDots]); colnames(groupingdotassignment) <- c(levelOfGroupingCaptial,levelOfDots)
      groupingdotassignment <- as.data.frame(groupingdotassignment)
      # if (anyDuplicated(colnames(groupingdotassignment)) != 0) {
      # }
      if (i == j) {
        colnames(groupingdotassignment)[2] <- paste0(colnames(groupingdotassignment)[1],"_too")
      }
      groupingdotassignment <- groupingdotassignment %>%
          group_by(!!sym(levelOfDots)) %>%
          summarise(new = first(!!sym(levelOfGroupingCaptial)))
      if (i == j) {
        colnames(groupingdotassignment)[2] <- paste0(colnames(groupingdotassignment)[1],"_too")
      } else if (i != j) {
              colnames(groupingdotassignment)[ncol(groupingdotassignment)] <- levelOfGroupingCaptial
      }
      
      # # then we sum all relative abundances within the current level of dots (e.g. genus), this will reorder the table alphabetically (not the samples!)
      ASV.by.sample.and.tax <- cbind(ASV.by.sample.and.tax[, 1:n.samples], ASV.by.sample.and.tax[, levelOfDots]); colnames(ASV.by.sample.and.tax)[ncol(ASV.by.sample.and.tax)] <- levelOfDots
      ASV.by.sample.and.tax <- ASV.by.sample.and.tax %>%
        group_by(!!sym(levelOfDots)) %>%
        summarise_all(.funs = sum)
      
      # # and then combine the two dataframes. this dataframe now has first all samples, then the level of dots (e.g. genus) and then the level of grouping (e.g. phylum)
      ASV.by.sample.and.tax <- cbind(ASV.by.sample.and.tax[,2:(n.samples+1)],groupingdotassignment)
      rm(groupingdotassignment)
    }
    
    # # 1. calcualte average relative abundance per clade (to scale points)
    asv.rowmeans <- rowMeans(ASV.by.sample.and.tax[,c(1:n.samples)])
    
    # # 2. calculate average of relative abundance per clade for each group (based on VisuaR Grouping1)
    averageGroup1 <- rowMeans(ASV.by.sample.and.tax[,which(M.groups==1)])
    averageGroup2 <- rowMeans(ASV.by.sample.and.tax[,which(M.groups==2)])
    averageGroup3 <- rowMeans(ASV.by.sample.and.tax[,which(M.groups==3)])
    rowmeans <- cbind(averageGroup1,averageGroup2,averageGroup3)
    
    # # 3. use result from 2. and divide by sum of relative abundances per clade over all groups (so they sum up to 1, rows sum up to 1 now) 
    normalized.rowmeans <- rowmeans/rowSums(rowmeans)
    rm(rowmeans,averageGroup1,averageGroup2,averageGroup3)
    colnames(normalized.rowmeans) <- M.projects.unique.ord
    
    # # 4. add average relative abundance per clade (from 1.)
    normalized.rowmeans <- cbind(normalized.rowmeans,asv.rowmeans); colnames(normalized.rowmeans)[ncol(normalized.rowmeans)] <- "mean"
    
    # # 5. extract names of top level of interest clades
    colored.clades <- levelOfGroupingTOP[c(1:numberofcoloredclades),1]
    
    # # change all ASVs not belonging to this to "Other"
    ASV.by.sample.and.tax <- ASV.by.sample.and.tax %>%
      mutate(Tax = ifelse(!!sym(levelOfGroupingCaptial) %in% colored.clades, !!sym(levelOfGroupingCaptial), "Other"))
    
    # # 6. create dataframe for ternary plot 
    ternary.data <- cbind(ASV.by.sample.and.tax$Tax,normalized.rowmeans); colnames(ternary.data)[1] <- "tax"; ternary.data <- as.data.frame(ternary.data)
    ternary.data[,2] <- as.numeric(ternary.data[,2])
    ternary.data[,3] <- as.numeric(ternary.data[,3])
    ternary.data[,4] <- as.numeric(ternary.data[,4])
    ternary.data[,5] <- as.numeric(ternary.data[,5])
    
    # # 7. create a color vector for the plot
    color_df <- data.frame(tax = unique(ternary.data$tax))
    if (any(color_df$tax == "Other")) {
      color_df$color <- ifelse(color_df$tax == "Other", "grey", colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(nrow(color_df) - 1))
      # color_df$color <- ifelse(color_df$tax == "Other", "grey", RColorBrewer::brewer.pal(n = nrow(color_df) - 1, name = "Dark2"))
    } else {
      if (nrow(color_df) == 2) {
        color_df$color <- RColorBrewer::brewer.pal(n = nrow(color_df), name = "Dark2")[1:2]
      } else {
        # color_df$color <- RColorBrewer::brewer.pal(n = nrow(color_df), name = "Dark2")
        color_df$color <- ifelse(color_df$tax == "Other", "grey", colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(nrow(color_df)))
        }
    }
    
    # # 8. Create the ternary plot with legends
    
    plot <- ggtern::ggtern(ternary.data, 
                           ggtern::aes(x=ternary.data[,2],
                                       y=ternary.data[,3],
                                       z=ternary.data[,4],
                                       size = mean, 
                                       color = tax)) +
      # ggplot2::scale_size_continuous(range = c(1, 7.5)) +  # Scale size by 0.5 times
      ggplot2::scale_color_manual(
        values = setNames(color_df$color, color_df$tax),
        name=levelOfGroupingCaptial)+
      # ggplot2::scale_color_brewer(palette = "Set1") +  # Use ColorBrewer palette
      ggplot2::labs(title=NULL,
                    x=colnames(ternary.data)[2],
                    y=colnames(ternary.data)[3],
                    z=colnames(ternary.data)[4],
                    size=paste0(levelOfDots," level\n","Relative Sequence\nAbundance"))+
      # ggtern::theme_rgbw() +
      ggtern::theme_void() +
      # ggtern::theme_custom(base_size = 10,
      #                      tern.plot.background = "transparent",
      #                      tern.panel.background = "transparent",
      #                      col.L = M.col[1],
      #                      col.T = M.col[2],
      #                      col.R = M.col[3],
      #                      col.grid.minor = NULL
      # ) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.direction = "vertical",  
                     legend.box = "vertical",      
                     legend.spacing.x = ggplot2::unit(0.2, "lines"),
                     legend.key.size = ggplot2::unit(0.001,"cm"))
    
    plot
    
    dir.create(file.path(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/TernaryPlots/")))
    pdf(file.path(paste0(VisuaRPath,"/",VisuaRName,"/Alpha_Diversity/TernaryPlots/",VisuaRName,"_",levelOfGrouping,"-",levelOfDots,"_ternary.pdf",sep="")),height=8,width=14,useDingbats=F); print(plot); dev.off()
    rm(plot)
  }

}
