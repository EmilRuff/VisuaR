![VisuaR Logo](https://user-images.githubusercontent.com/62703732/82889437-41dacf80-9f3a-11ea-8671-8296d87e47d2.png)

# 
R-based workflow for the analysis and visualization of marker gene data

VisuaR's main goal is to provide a set of analyses and visualizations that are often used in the analyses of marker gene data. Although, the script was originally designed to analyze 16S rRNA gene amplicon sequence data processed by dada2 it can accomodate marker gene data from diverse sources (e.g., metagenomics, long-read sequencing) as long as the input data is correctly formatted. The VisuaR workflow is based on commands from other R packages (e.g. vegan, stringr, plyr, ggplot2), as well as on custom-build scripts. It takes dada2 files (https://github.com/benjjneb/dada2) as input and generates (near-)publication ready figures as output. In addition to NMDS ordinations, cluster dendrograms and UpSet diagrams the workflow features custom scripts that calculate the percentage of singletons or shared ASVs, and visualize diversity indices, or the relative sequence abundance of clades at specific taxonomic levels. All analyses and visualizations can be performed using a contextual data matrix to create figures colored by groupings (e.g. environmental parameters).

## Input Files
1.  A sequence table without chimeras as produced by dada2 after using dada2::makeSequenceTable and dada2:: removeBimeraDenovo (Samples by ASVs)
2.  A taxonomy file as produced by dada2:: assignTaxonomy or dada2::addSpecies (ASVs by Taxonomy)
3.	Optional: A contextual data file (.txt)
        a.	First column: sample names (the same names as used for the dada2 analysis output. Those names are not necessarily the same names as used during the sequencing as they might change during the dada2 analysis depending on your input)
        b.	Other columns containing your contextual data. 
        c.	Replace missing values by NA (no empty cells!)
        d.	Avoid spaces and special symbols, replace them by underscores if necessary

## Output Files
1. Alpha Diversity:
        • Composition Plots on all taxonomic levels (Bubble and Box)
        • Plots of Diversity Indices (Box and Violin)
        • Species Accumulation Curves
        • Rarefaction Curve 
        • Composition Tables on all taxonomic levels
        • Relative abundance, abundance, presence absence
        • Diversity Indices table
        • ASV by Sample table with abundances
        • ASV by sample table with presences/absences
        • ASV by Sample table with relative abundances
        • Reads and ASVs per sample table
        • Reads per ASV table
2. Beta Diversity:
        • NMDS plots
        • Anosim plot
        • UpsetR plot
        • VennDiagram
        • Dendrograms (not optimized)
        • .txt files with
        • Percent shared ASVs between groups
        • Shared ASVs between groups
        • Total ASVs between groups
        • Unique ASVs between groups
        • Unique ASVs per grouping
3. Analysis Log
4. Session Info (R version and loaded packages)
5. R workspace
6. ASV summary (a list of the ASV identifiers, the assigned numbers in VisuaR and the respective taxonomy)
