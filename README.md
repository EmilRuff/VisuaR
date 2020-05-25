![VisuaR Logo](https://user-images.githubusercontent.com/62703732/82836047-4f9a4180-9eb5-11ea-987d-b67a5b56241e.png)

# VisuaR
R-based workflow for the analysis and visualization of dada2-processed sequencing data

VisuaR's main goal is to provide a set of analyses and visualizations that are often used in the analyses of amplicon sequencing data. The VisuaR workflow is based on commands from other R packages (e.g. vegan, stringr, plyr, ggplot2), as well as on custom-build scripts. It takes dada2 files (https://github.com/benjjneb/dada2) as input and generates (near-)publication ready figures as output. In addition to NMDS ordinations, cluster dendrograms and UpSet diagrams the workflow features custom scripts that calculate the percentage of singletons or shared ASVs, and visualize diversity indices, or the relative sequence abundance of clades at specific taxonomic levels. All analyses and visualizations can be performed using a contextual data matrix to create figures colored by groupings (e.g. environmental parameters).

## Input Files
1.  A sequence table without chimeras as produced by dada2 after using dada2::makeSequenceTable and dada2:: removeBimeraDenovo (Samples by     ASVs)
2.  A taxonomy file as produced by dada2:: assignTaxonomy or dada2::addSpecies (ASVs by Taxonomy)
3.	Optional: A contextual data file (.txt)
        a.	First column: sample names (the same names as used for the dada2 analysis output. Those names are not necessarily the same                 names as used during the sequencing as they might change during the dada2 analysis depending on your input)
        b.	Other columns containing your contextual data. 
        c.	Replace missing values by NA (no empty cells!)
        d.	Avoid spaces and special symbols, replace them by underscores if necessary
