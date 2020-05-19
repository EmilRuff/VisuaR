# #===== Calculate Procrustes ordination ================================
# # find a really good introduction into Procrustes and PROTEST here: http://jackson.eeb.utoronto.ca/procrustes-analysis/
# # RDocumentation: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/procrustes
# # Procrustes (least-squares orthogonal mapping) compares two sets of configurations (e.g. NMDS ordinations)
# # Tries to minimize the 'sum of squared deviations' (=error, m^2) through translating, rotating and dilating one configuration to match the other configuration (target)
# # m^2 is a measure of fit but does not include information on statistical significance of the solution (is the solution better than a solution by chance). Significance testing is performed by PROTEST
# # The deviations between two data points are called residuals. Small residual = close agreement between corresponding data points
# # This is useful to evaluate e.g. how Singletons or contaminants influence the community structure.
# # High congruence indicates high community similarity.

# # PROTEST: PROcrustean randomization TEST
# # Permutation approach to test whether m^2 is smaller than expected by chance
# # compares whether the original m^2 is smaller than or equal to the m^2 value obtained from the fit of the randomized dataset to the second data set (permustats())
# # it tabulates the number of times where the observed m^2 value was smaller than or equal to that obatained from the randomized dataset
library(vegan)

NameForProcrustesAnalysis='Name'
FolderForProcrustesPlot=file.path('')                           # The directory where you want to save your plots

# # Input from corresponding VisuaR analysis
M.colvec=readRDS('C:/.../Mcolvec.rds')                          # this .rds file can be found in your project folder ('VisuarProjectName_Mcolvec.rds')
M.projects.unique.ord=readRDS('C:/.../Mprojectsuniqueord.rds')  # this .rds file can be found in your project folder ('VisuarProjectName_Mprojectsuniqueord.rds')
M.col=readRDS('C:/.../Mcol.rds')                                # this .rds file can be found in your project folder ('VisuarProjectName_Mcol.rds')
TargetMatrix=readRDS("C:/.../ProjectName_distanceMatrix.rds")   # The target matrix is the original matrix, the 'MatrixtobeRotated' will be scaled and rotated to fit to the target matrix. This matrix can be found in your analysis folder under Beta_Diversity, 'ProjectName_distanceMatrix.rds'
MatrixtoBeRotated=readRDS("C:/.../ProjectName2_distanceMatrix.rds")

# # PROTEST analysis is a modified version of procrustes that has more versatility and statistical testing
M.prot=protest(TargetMatrix,MatrixtoBeRotated,permutations=1000) # tests the non-randomness (significance) between two configurations

pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes_PROTEST.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
plot(M.prot,ylab = '',xlab = '',ar.col = '#595959',to.target = F,main = 'Procrustes Analysis (PROTEST)')  
points(M.prot,col=M.colvec,pch=19,display = 'target')
par(new=T,xpd=T)                                                                                                       
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                                            
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.46,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.31,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.21,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.12,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()

pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes_PROTEST_names.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
plot(M.prot,ylab = '',xlab = '',ar.col = '#595959',to.target = F,main = 'Procrustes Analysis (PROTEST)')  
par(mar=c(2, 2, 2, 7), xpd=F)                                                                           
text(M.prot,col=M.colvec,display = 'target',cex=0.5)
par(new=T,xpd=T)                                                                                                       
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                                            
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.46,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.31,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.21,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.12,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()

pdf(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_Procrustes_PROTEST_residuals.pdf',sep='')),height=5,width=6,useDingbats=F)
par(mar=c(7, 4, 2, 7), xpd=F)                                                                   # mar: margin sizes in the following order: bottom, left, top, and right. There is more space at the right hand site for the figure legend, xpd='T': all plotting is clipped to the figure region (not only to the plot region). This allows to place figure legens outside of the plot
plot(M.prot,main = 'Procrustes Analysis - Residuals',kind = 2,xaxt='n',xlab='',ylab = 'Residuals',col=M.colvec,lwd=3)  # plot ordination
axis(1, at=1:length(M.colvec),labels =rownames(TargetMatrix$points),las=3)
par(new=T,xpd=T)                                                                                # only using par(new=T) it is possible to plot the figure legend independent of the axis scaling to a fixed place inside the pdf.
plot.new()
legend(x=1.01,y=1.05,legend='Groups',bty='n',cex=0.75,text.font = 2)                            # text.font=2: prints text in bold, bty='n': no box will be drawn around the legend. Alternative: bty='o'
legend(1.045,1.00,legend=c(M.projects.unique.ord),fill=c(M.col),bty='n',cex=0.75)
legend(1.01,0.46,legend=paste0("Procrustes\nSum of squared\ndeviations:"),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.27,legend=paste0(round(M.prot$ss,3)),bty='n',cex=0.75)
legend(1.01,0.21,legend=paste0('Protest\nSignificance'),bty='n',cex=0.75,text.font = 2)
legend(1.01,0.09,legend=paste0(round(M.prot$signif,3)),bty='n',cex=0.75)
dev.off()

permustats(M.prot)
summary(M.prot)
residuals(M.prot)

save.image(file.path(FolderForProcrustesPlot,paste(NameForProcrustesAnalysis,'_procrustes_workspace','.RData',sep='')))
