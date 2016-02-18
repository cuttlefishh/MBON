# MBON_methods_experiment
R script to process OTU_table file output by banzai pipeline. This script will be compatible for all target loci (e.g. 12S, 16S, 18S, 28S). Currently the script includes the following input files:

1. OTU_table
2. Metadata file
3. MEGAN csv file containing 'No hits'

And following steps:

1. Removal of contaminant reads identified by positive and/or negative controls
2. Removal of OTUs assigned as 'no hits' via BLAST
3. Data normalization with DESeq2
4. Data analyses and plots
    a. Species richness
    b. Multiple linear regression
    c. NMDS/hierarchical clustering

Analyses are at the OTU level but we can work to incorporate taxonomy as well. I am also working to incorporate species occupancy modeling based on script provided by Lahoz-Monfort et al. 2015 and Ryan Kelly. This will identify PCR/sequencing errors (i.e. low frequency noise), and take the place of ad hoc subtraction of reads/OTUs/taxon that occur in X out of Y replicates.
