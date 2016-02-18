#########################################################################################################
############################################################################################################################################
#This script analyzes data from the MBON extraction/filtration methods experiment
#It compares species richness and abundance across different extraction methods and filter types
############################################################################################################################################
############################################################################################################################################
rm(list = ls())
library(DESeq2)
library(vegan)
library(ggplot2)

OTU_table="/Users/jport/Desktop/Analysis_20160202_2050_SERVER/OTU_table.csv" #12S
#OTU_table="/Users/jport/Desktop/18S_kelp_MBON_methods/OTU_table.csv" #18S
metadata="/Users/jport/Desktop/MiSeq_MBON_ExtFilt_MontBay_Metadata.csv" #12S
#metadata="/Users/jport/Desktop/18S_kelp_MBON_methods/BOG_kelp18S_metadata.csv" #18S
megan_nohits_file="/Users/jport/Desktop/Analysis_20160202_2050_SERVER/Took_OTUs.fasta_from_server_and_BLAST_separately_on_server/megan_nohits.csv" #12S
#megan_nohits_file="/Users/jport/Desktop/Analysis_20160202_2050_SERVER/Took_OTUs.fasta_from_server_and_BLAST_separately_on_server/megan_nohits.csv" #18S
metadata_col = 12 #12S
#metadata_col = 16 #18S

##########################################################################################
# 1. REMOVE CONTAMINATION
##########################################################################################
OTUs = read.csv(OTU_table, header=TRUE, row.names=1)
meta=read.csv(metadata, header=TRUE)

env = meta[meta$sample_type=="Environmental",]
env_tags = paste0("lib_", env$library,"_tag_", env$tag_sequence)
samples_env <- OTUs[,match(env_tags,names(OTUs))]

NEGCONTROL=meta[meta$sample_type=="NegControl",]  #select relevant samples from the metadata
NEGCONTROL_tags=paste0("lib_", NEGCONTROL$library,"_tag_", NEGCONTROL$tag_sequence)  #generate tag-ids for samples from relevant experiment
samples_negcontrol <- OTUs[,match(NEGCONTROL_tags,names(OTUs))]
samples_negcontrol= samples_negcontrol[colSums(samples_negcontrol)>10000] #only use the NTC samples with a reasonable number of reads sequenced.

POSCONTROL=meta[meta$sample_name=="Swordfish",]
POSCONTROL_tags=paste0("lib_", POSCONTROL$library,"_tag_", POSCONTROL$tag_sequence)  #generate tag-ids for samples from relevant experiment
samples_poscontrol <- OTUs[,match(POSCONTROL_tags,names(OTUs))]

# IDENTIFY CONTROL TAXON (Only if including a positive control; need to modify code if using artificial community instead of tissue sample)
# Which OTU is (must be a row name in the samples_all matrix)
# control_taxon <- "OTU_1"
control_taxon <- names(which.max(rowSums(samples_poscontrol)))

contam_proportion_df=data.frame(samples_negcontrol) #Include positice control here too as in Ryan's code?
contam_proportion_df=scale(contam_proportion_df, center=F, scale=colSums(contam_proportion_df))  #for each dup/OTU in controls, calc frequency
contam_proportion_vec =apply(contam_proportion_df, 1, max)  #take max of these frequencies.  this is the max proportional contribution of contamination to each of these OTUs/dups in the dataset. Accordingly, you will then subtract this proportion from each dup/OTU occurrence in the dataset.  A different (better?) way to do this would be to fit a likelihood distribution to the observed proportions of the dups/OTUs in the controls, and come up with the most likely proportion attributable to contamination.  But, in most cases, the most likely contribution is at or near zero (using a beta distrib), and in many cases the models don't converge anyway.

#Step 1: remove the fraction of each OTU abundace attributable to contamination
DECONTAM = samples_env
for (i in 1:ncol(DECONTAM)){
  nreads=sum(DECONTAM[,i])
  contam_proportion_vec_forSample= ceiling(contam_proportion_vec*nreads)
  DECONTAM[,i]= samples_env[,i]-contam_proportion_vec_forSample
}
DECONTAM[DECONTAM <0]<-0  #assign zero to any reads for which there are a negative (corrected) number of reads
DECONTAM = DECONTAM[rowSums(DECONTAM)>0,]  #remove reads that no longer occur in the data

#Step 2: remove control taxon OTU from environmental samples in which it occurs #Only if including a positive control
DECONTAM <- DECONTAM[!(rownames(DECONTAM) %in% control_taxon),]

###############################################################################################
# 2. REMOVE OTUs assigned to "No hits" (If analyzing data at OTU level)
###############################################################################################
OTUs2 = DECONTAM
#OTUs = read.csv("/Users/jport/Desktop/Analysis_20160202_2050_SERVER/OTU_table.csv", header=TRUE, row.names=1)
megan_nohits = read.csv(megan_nohits_file, header=FALSE)
OTUs2 = subset(OTUs2, !(row.names(OTUs2) %in% megan_nohits[,1]))

#################################################################################################################
# 3. NORMALIZE DATASET WITH DESEQ2
#################################################################################################################
#Load count datasheet
deseq_count=OTUs2 #(rows=samples, columns=taxa, cells=sequence read counts)

#Remove any samples whose total read count = 0 (incompatible with DESeq2 normalization)
deseq_count=OTUs2[colSums(OTUs2) != 0] #remove samples whose total counts = 0 (incompatible with size factor estimation)

#Load metadata sheet
deseq_col=env
#deseq_col=read.csv("/Users/jport/Desktop/MiSeq_MBON_ExtFilt_MontBay_Metadata.csv", sep=",", header=TRUE)

          #If any samples had total taxon count  = 0 per above, subtract these samples from metadata sheet (required to match number of samples in count datasheet) (need to make this code more flexible, using paste lib tag?)
                #e.g. deseq_col=deseq_col[(deseq_col$Library != "lib1" | deseq_col$TAG_SEQUENCE != "AGCCTC"),]

#Create DESeqDataSet object
deseq_obj=DESeqDataSetFromMatrix(deseq_count, deseq_col, design= ~ sample_name)

#Estimate size factors when zeros are present in each row
fordeseq=counts(deseq_obj)
geoMeans = apply(fordeseq, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
deseq_obj= estimateSizeFactors(deseq_obj, geoMeans=geoMeans)

#Format normalized count table
deseqtable=counts(deseq_obj, normalized=TRUE) #store normalized counts table
    colnames(deseqtable)=colnames(deseq_count) #add back column names to deseq table
    deseqtable=t(deseqtable) #transpose deseq table so rows = samples and columns = taxa
    deseqtable=as.data.frame(deseqtable)

#################################################################################################################
# 4. TRANSFORM COUNT DATA TO BINARY AND CALCULATE SPECIES RICHNESS
#################################################################################################################
datasheet_binary = deseqtable
datasheet_binary[datasheet_binary>0] <-1
#datasheet_binary <-specnumber(datasheet_binary) #vegan code for species richness

#sum binary count by sample (i.e. row)
binsum=apply(datasheet_binary, 1, sum)

#Append row sums to binary table
datasheet_binary = cbind(datasheet_binary, binsum)

datasheet_binary = cbind(env, datasheet_binary) #make sure sample order in datasheet = sample order in metadata sheet
#datasheet_binary = subset(datasheet_binary, sample_type=="Environmental") #remove positive and negative controls

#################################################################################################################
#5. GENERATE PLOTS COMPARING SPECIES RICHNESS AND ASSOCIATED STATS
#################################################################################################################
#by extraction method
ggplot(datasheet_binary, aes(x=extraction_method, y=binsum))+geom_boxplot()+theme_classic()+ylab("No. of OTUs")+xlab("Extraction method")

        #generate means and standard deviations by extraction method
        extraction_means=aggregate(datasheet_binary$binsum~datasheet_binary$extraction_method, FUN=mean)
        extraction_sd=aggregate(datasheet_binary$binsum~datasheet_binary$extraction_method, FUN=sd)
        
        #statistical test: does species richness differ between extraction methods?
        kruskal.test(binsum~extraction_method, data=datasheet_binary) #is sample size large enough to run this test?
        
#by filter type
ggplot(datasheet_binary, aes(x=filter_type, y=binsum))+geom_boxplot()+theme_classic()+ylab("No. of OTUs")+xlab("Filter type")

        #generate means and standard deviations by filter type
        filter_means=aggregate(datasheet_binary$binsum~datasheet_binary$filter_type, FUN=mean)
        filter_sd=aggregate(datasheet_binary$binsum~datasheet_binary$filter_type, FUN=sd)

        #statistical test: does species richness differ between filter types? 
        kruskal.test(binsum~filter_type, data=datasheet_binary) #is sample size large enough to run this test?
        
#by filter type within extraction method
ggplot(datasheet_binary, aes(x=extraction_method, y=binsum, fill=filter_type))+geom_boxplot()+theme_classic()+ylab("No. of OTUs")+xlab("Extraction method")+theme(legend.position=c(.88,.8))+theme(legend.text=element_text(size=10))+theme(legend.key.size = unit(.5, "cm"))

        #Statistical tests to compare methods (mann-whitney to compare 2 groups, kruskal wallis to compare 3 or more groups) 
        kruskal.test(dneasy$binsum~filter_type, data=dneasy) #by filter type within extraction method 
        kruskal.test(binsum~filter_type, data=mobio) #by filter type within extraction method 
        kruskal.test(binsum~filter_type, data=pc) #by filter type within extraction method

#by dilution factor
boxplot(binsum~dilution, data=datasheet_binary)

        #generate means and standard deviations by dilution factor
        dilution_means = aggregate(binsum~dilution, data=datasheet_binary, FUN=mean)
        dilution_sd = aggregate(binsum~dilution, data=datasheet_binary, FUN=mean)

        #statistical test: does species richness differ between dilution factors? 
        kruskal.test(binsum~dilution, data=datasheet_binary) #is sample size large enough to run this test?
        
#Hierarchical clustering
dist=vegdist(datasheet_binary[,c(metadata_col:ncol(datasheet_binary))], method="bray")
clust.res=hclust(dist, method="average")
plot(clust.res, labels=datasheet_binary$extraction_method, hang=-1, cex=0.5, main=NULL, xlab="Sample") #or use labels=n if this doesn't work correctly

#NMDS
dna.mds <- metaMDS(datasheet_binary[,c(metadata_col:ncol(datasheet_binary))], distance="bray", trace = FALSE, trymax=100) #Multi-dimensional Scaling (MDS)
    
    ordiplot(dna.mds,type="n", xlim=c(-1,1))

          for (i in 1:length(unique(datasheet_binary$extraction_method))){
              ordihull(dna.mds, groups=datasheet_binary$extraction_method, draw="lines", show.groups = unique(datasheet_binary$extraction_method)[i] ,col=c(rainbow(8, alpha = 1))[i],label=F)
          }

          orditorp(dna.mds,display="sites")
          

#Subset data table by extraction method
dneasy=subset(datasheet_binary, datasheet_binary$extraction_method=="DNEASY")
mobio=subset(datasheet_binary, datasheet_binary$extraction_method=="MOBIO")
pc=subset(datasheet_binary, datasheet_binary$extraction_method=="PC")




#################################################################################################################
#GENERATE PLOTS AND STATS COMPARING COMMUNITY COMPOSITION ACROSS METHODS
#################################################################################################################

datasheet=cbind(sample_metadata[,1:8], deseqtable) 
adonis(datasheet[,c(9:ncol(datasheet))] ~ datasheet$extraction_method, method="bray", perm=500) #make sure no NA in datasheet

#################################################################################################################
#UNUSED CODE
#################################################################################################################

#make table of only taxon counts
taxontable_counts=taxontable[,c(9:length(taxontable))]

#plot distribution of counts by sample type (e.g. env, pos, neg)
boxplot(sumrows~Sample_Type, data=taxontable)

#Plot total DNA concentrations
##grouped by extraction method
ggplot(na.omit(dna2), aes(x = ext, y= conc, fill=location)) + geom_boxplot()+scale_fill_manual(values = rep(k, 3))+xlab("Extraction Method")+ylab("DNA Concentration (ng/ul)")+theme_bw()+theme(legend.position = c(0.8, .91))+theme(axis.text=element_text(size=17), legend.text=element_text(size=18), text=element_text(size=16))
##grouped by location
ggplot(na.omit(dna2), aes(x = location, y= conc, fill=ext)) + geom_boxplot()+scale_fill_manual(values = rep(k,2))+xlab("Sample")+ylab("Total DNA Concentration (ng/ul)")+theme_bw()+theme(legend.position = c(0.8, .85))+theme(axis.text=element_text(size=18), legend.text=element_text(size=20), text=element_text(size=20), axis.title.x=element_text(vjust=0.0001), axis.title.y=element_text(vjust=1.3))