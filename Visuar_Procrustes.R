# #===== Calculate Procrustes ordination ================================
# # RDocumentation: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/procrustes
# # find a really good introduction into Procrustes and PROTEST here: http://jackson.eeb.utoronto.ca/procrustes-analysis/

# # Procrustes (least-squares orthogonal mapping): rotates configuration to maximum similarity with another configuration
# # typically used in comparison of ordination results. It is particularly useful in comparing alternative solutions in multidimensional scaling
# # # Tries to minimize the 'sum of squared deviations' (=error, m^2) through translating, rotating and dilating one configuration to match the other configuration (target)
# # # m^2 is 0 for same matrices 
# # # # m^2 is a measure of fit but does not include information on statistical significance of the solution
# # # # A smaller SSD (sum of squared deviations) indicates a better alignment between the two sets of points

# # # Residuals
# # # # The deviations between two data points are called residuals. Small residual = close agreement between corresponding data points
# # # # residuals let us identify the samples with the worst fit. horizontal lines are 25, 50 and 75 % quanitles of the residuals

# # Protest
# # Procrustean randomization test: tests the non-randomness (significance) between two configurations
# # # Permutation approach to test whether m^2 is smaller than expected by chance
# # # compares whether the original m^2 is smaller than or equal to the m^2 value obtained from the fit of the randomized dataset to the second data set (permustats())
# # # it tabulates the number of times where the observed m^2 value was smaller than or equal to that obatained from the randomized dataset
# # # it is a permutational test of the significance of symmetric procrustes analyses (symmetric: order of origin and target matrix do not matter)
# # # be aware that as this is a permutation test, the significance can not be less than 1/(permutations + 1)

# # add info
# # metaMDS() is also doing procrustes
# # # two final configurations are considered to have converged (arrived at essentially the same solution) when the rmse (root mean squared error) is less than 0.01, and no single residual value exceeds 0.005. (https://www.flutterbys.com.au/stats/tut/tut15.1.html)

# # Input
# # we are using the _distanceMatrix2d.rds from VisuaR which is the results of M.mMDS=metaMDS(M.dist)


library(vegan)

NameForProcrustesAnalysis=''

# # Input from corresponding VisuaR analysis #####
VisuaRname.origin <- "" # VisuaR name of the "origin" distance matrix(the target will be scaled and rotated to fit on the original matrix)
path.to.origin <- file.path("") 

VisuaRname.target <- ""
path.to.target <- file.path("")

FolderForProcrustesPlot=file.path("") # The directory where you want to save your plots

dir.create(FolderForProcrustesPlot)
# # read in of files ####

M.colvec <- readRDS(file.path(path.to.origin,VisuaRname.origin,paste(VisuaRname.origin,'_Mcolvec.rds',sep=''))) # this .rds file can be found in your project folder ('VisuarProjectName_Mcolvec.rds')
M.projects.unique.ord <- readRDS(file.path(path.to.origin,VisuaRname.origin,paste(VisuaRname.origin,'_Mprojectsuniqueord.rds',sep=''))) # this .rds file can be found in your project folder ('VisuarProjectName_Mprojectsuniqueord.rds')
M.col <- readRDS(file.path(path.to.origin,VisuaRname.origin,paste(VisuaRname.origin,'_Mcol.rds',sep=''))) # this .rds file can be found in your project folder ('VisuarProjectName_Mcol.rds')
TargetMatrix <- readRDS(file.path(path.to.origin,VisuaRname.origin,"Beta_Diversity",paste(VisuaRname.origin,'_distanceMatrix2d.rds',sep='')))

MatrixtoBeRotated <- readRDS(file.path(path.to.target,VisuaRname.target,"Beta_Diversity",paste(VisuaRname.target,'_distanceMatrix2d.rds',sep='')))

# # PROTEST analysis ####
M.prot=protest(TargetMatrix,MatrixtoBeRotated,permutations=1000) 
M.prot 


# # Plots ####
pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
plot(M.prot,ylab = '',xlab = '',ar.col = '#595959',to.target = F,main = paste0("Procrustes - ",VisuaRname.origin," vs ", VisuaRname.target),cex.main=0.6)  
points(M.prot,col=M.colvec,pch=19,display = 'target')
par(new=T,xpd=T)                                                                                                       
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                                            
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.64,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.50,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.46,legend=paste0("Root mean\nsquared errors\n(rmse):"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.32,legend=paste0(round(as.numeric(summary(M.prot) [6]),3)),bty='n',cex=0.75)
legend(1.01,0.26,legend=paste0("Residuals\nmin,mean,max:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.18,legend=paste0(round(min(residuals(M.prot)),3),", ",round(mean(residuals(M.prot)),3),", ",round(max(residuals(M.prot)),3)),bty='n',cex=0.75)
legend(1.01,0.12,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.04,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()


pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes_names.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
plot(M.prot,ylab = '',xlab = '',ar.col = '#595959',to.target = F,main = paste0("Procrustes - ",VisuaRname.origin," vs ", VisuaRname.target),cex.main=0.6)  
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
text(M.prot,col=M.colvec,display = 'target',cex=0.5)
par(new=T,xpd=T)                                                                                                       
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                                            
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.64,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.50,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.46,legend=paste0("Root mean\nsquared errors\n(rmse):"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.32,legend=paste0(round(as.numeric(summary(M.prot) [6]),3)),bty='n',cex=0.75)
legend(1.01,0.26,legend=paste0("Residuals\nmin,mean,max:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.18,legend=paste0(round(min(residuals(M.prot)),3),", ",round(mean(residuals(M.prot)),3),", ",round(max(residuals(M.prot)),3)),bty='n',cex=0.75)
legend(1.01,0.12,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.04,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()

min(residuals(M.prot))

pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes_residuals.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(3, 4, 2, 7), xpd=F)                                                                   # mar: margin sizes in the following order: bottom, left, top, and right. There is more space at the right hand site for the figure legend, xpd='T': all plotting is clipped to the figure region (not only to the plot region). This allows to place figure legens outside of the plot
plot(M.prot,main = paste0("Procrustes residuals - ",VisuaRname.origin," vs ", VisuaRname.target),kind = 2,xaxt='n',xlab='',ylab = 'Residuals',col=M.colvec,lwd=3, cex.main= 0.6)  # plot ordination
axis(1, at=1:length(M.colvec),labels =rownames(TargetMatrix$points),las=3,lwd.ticks = 0.2,cex.axis=0.1)
par(new=T,xpd=T)                                                                                # only using par(new=T) it is possible to plot the figure legend independent of the axis scaling to a fixed place inside the pdf.
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                            # text.font=2: prints text in bold, bty='n': no box will be drawn around the legend. Alternative: bty='o'
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.68,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.53,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.49,legend=paste0("Root mean\nsquared errors\n(rmse):"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.34,legend=paste0(round(as.numeric(summary(M.prot) [6]),3)),bty='n',cex=0.75)
legend(1.01,0.28,legend=paste0("Residuals\nmin,mean,max:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.19,legend=paste0(round(min(residuals(M.prot)),3),", ",round(mean(residuals(M.prot)),3),", ",round(max(residuals(M.prot)),3)),bty='n',cex=0.75)
legend(1.01,0.13,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.04,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()

sink(file=(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_log.txt',sep=''))),append=TRUE,type='output')

cat('Procrustes Analysis:',"\n",sep='')
cat(file.path(FolderForProcrustesPlot,NameForProcrustesAnalysis),"\n\n", sep="")

cat("Target matrix and info from VisuaR analysis: ","\n", sep="")
cat(file.path(path.to.origin,VisuaRname.origin),"\n\n", sep="")

cat("Matrix to be rotated to fit the target matrix VisuaR analysis: ","\n", sep="")
cat(file.path(path.to.target,VisuaRname.target),"\n\n", sep="")

cat("Results of Protest analysis: ","\n", sep="")
M.prot
permustats(M.prot)
summary(M.prot)

sink()

data.table::fwrite(as.data.frame(residuals(M.prot)),file=file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_residuals.csv',sep='')),col.names = T,row.names = T)
# residuals(M.prot)

# save.image(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_procrustes_workspace','.RData',sep='')))

closeAllConnections() # closes all currently open connections.

rm(list=ls())
gc()
