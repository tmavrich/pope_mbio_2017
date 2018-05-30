#Script to analyze Actinobacteriophage_789 and Cyanobacteriophage_209 gene content data
#as published in Pope et al., mBio, 2017.



###Define functions

#Plot gene content dissimilarity values by phage
plot_gcd_values_by_phage <- function(input_table,phage_name){
  
  reduced_table <- subset(input_table,input_table$reference == phage_name | input_table$query == phage_name)
  reduced_table_sorted <- reduced_table[order(reduced_table$pham_pham_dissimilarity,decreasing=TRUE),]
  
  par(mar=c(4,8,4,4))
  plot(reduced_table_sorted$pham_pham_dissimilarity,ylim=c(0,1),pch=16,cex=3,cex.axis=2,ann=FALSE,main=NULL,las=1,frame.plot = FALSE)

  
}





###Import Cyano209 data

#Import list of unique (non-redundant) pairwise comparisons
cyano209_unique_pairwise_comparisons <- read.csv("cyano209_pairwise_comparison_list.csv",sep=",",header=TRUE)


#Import host taxonomy data
cyano209_metadata_table <- read.csv("cyano209_phage_metadata.csv",sep=",",header=TRUE)
cyano209_metadata_table[cyano209_metadata_table == "Unspecified"] <- NA

#Merge with comparison list
names(cyano209_metadata_table) <- c("query_phagename","query_host_genus","query_genome_size","query_accession","query_gregory2016","query_site","query_lineage")
cyano209_main_data_table <- merge(cyano209_unique_pairwise_comparisons,cyano209_metadata_table,by.x="query",by.y="query_phagename")
names(cyano209_metadata_table) <- c("ref_phagename","ref_host_genus","ref_genome_size","ref_accession","ref_gregory2016","ref_site","ref_lineage")
cyano209_main_data_table <- merge(cyano209_main_data_table,cyano209_metadata_table,by.x="reference",by.y="ref_phagename")
names(cyano209_metadata_table) <- c("phagename","host_genus","genome_size","accession","gregory2016","site","lineage")



#Import pham data and merge with comparison list
cyano209_pham_table <- read.csv("cyano209_pairwise_pham_proportions.csv",sep=",",header=TRUE)
cyano209_pham_table$pham_dissimilarity <- 1 - cyano209_pham_table$average_shared_proportion

names(cyano209_pham_table) <- c("pham_phage1","pham_phage1_number_unshared_phams","pham_phage1_shared_proportion",
                             "pham_phage2","pham_phage2_number_unshared_phams","pham_phage2_shared_proportion",
                             "pham_number_shared_phams","pham_average_shared_proportion",
                             "pham_jaccard_similarity","pham_shared_pham_distribution_mean",
                             "pham_shared_pham_distribution_median","pham_shared_pham_distribution_max",
                             "pham_unshared_pham_distribution_mean","pham_unshared_pham_distribution_median",
                             "pham_unshared_pham_distribution_max","pham_unshared_orpham_count",
                             "pham_pham_dissimilarity")
                             
cyano209_pham_table$pham_phage1_phage2 <- paste(cyano209_pham_table$pham_phage1,"_",cyano209_pham_table$pham_phage2,sep="")
cyano209_pham_table$pham_phage1_phage2 <- as.factor(cyano209_pham_table$pham_phage1_phage2)
cyano209_main_data_table <- merge(cyano209_main_data_table,cyano209_pham_table,by.x="ref_query",by.y="pham_phage1_phage2")








#Compare ref and query host data columns
cyano209_main_data_table$host_genus_compare <- ifelse(cyano209_main_data_table$ref_host_genus==cyano209_main_data_table$query_host_genus,
                                                      as.character(cyano209_main_data_table$ref_host_genus),"different")

cyano209_main_data_table$gregory2016_compare <- ifelse(cyano209_main_data_table$ref_gregory2016==cyano209_main_data_table$query_gregory2016,
                                                       as.character(cyano209_main_data_table$ref_gregory2016),"different")

#Convert to factor class
cyano209_main_data_table$host_genus_compare <- as.factor(cyano209_main_data_table$host_genus_compare)
cyano209_main_data_table$gregory2016_compare <- as.factor(cyano209_main_data_table$gregory2016_compare)

#Convert all Unspecified fields to NA missing value
cyano209_main_data_table[cyano209_main_data_table == "Unspecified"] <- NA






###Import Actino789 data

#Import list of unique (non-redundant) pairwise comparisons
actino789_unique_pairwise_comparisons <- read.csv("actino789_pairwise_comparison_list.csv",sep=",",header=TRUE)


#Import host taxonomy data
actino789_metadata_table <- read.csv("actino789_phage_metadata.csv",sep=",",header=TRUE)
actino789_metadata_table[actino789_metadata_table == "Unspecified"] <- NA

#Merge with comparison list
names(actino789_metadata_table) <- c("query_phagename","query_host_genus","query_cluster","query_subcluster","query_genome_size")
actino789_main_data_table <- merge(actino789_unique_pairwise_comparisons,actino789_metadata_table,by.x="query",by.y="query_phagename")
names(actino789_metadata_table) <- c("ref_phagename","ref_host_genus","ref_cluster","ref_subcluster","ref_genome_size")
actino789_main_data_table <- merge(actino789_main_data_table,actino789_metadata_table,by.x="reference",by.y="ref_phagename")
names(actino789_metadata_table) <- c("phagename","host_genus","cluster","subcluster","genome_size")




#Import pham data and merge with comparison list
actino789_pham_table <- read.csv("actino789_pairwise_pham_proportions.csv",sep=",",header=TRUE)
actino789_pham_table$pham_dissimilarity <- 1 - actino789_pham_table$average_shared_proportion

names(actino789_pham_table) <- c("pham_phage1","pham_phage1_number_unshared_phams","pham_phage1_shared_proportion",
                                "pham_phage2","pham_phage2_number_unshared_phams","pham_phage2_shared_proportion",
                                "pham_number_shared_phams","pham_average_shared_proportion",
                                "pham_jaccard_similarity","pham_shared_pham_distribution_mean",
                                "pham_shared_pham_distribution_median","pham_shared_pham_distribution_max",
                                "pham_unshared_pham_distribution_mean","pham_unshared_pham_distribution_median",
                                "pham_unshared_pham_distribution_max","pham_unshared_orpham_count",
                                "pham_pham_dissimilarity")



actino789_pham_table$pham_phage1_phage2 <- paste(actino789_pham_table$pham_phage1,"_",actino789_pham_table$pham_phage2,sep="")
actino789_pham_table$pham_phage1_phage2 <- as.factor(actino789_pham_table$pham_phage1_phage2)
actino789_main_data_table <- merge(actino789_main_data_table,actino789_pham_table,by.x="ref_query",by.y="pham_phage1_phage2")





#Compare ref and query host data columns
actino789_main_data_table$host_genus_compare <- ifelse(actino789_main_data_table$ref_host_genus==actino789_main_data_table$query_host_genus,
                                                       as.character(actino789_main_data_table$ref_host_genus),"different")

actino789_main_data_table$cluster_compare <- ifelse(actino789_main_data_table$ref_cluster==actino789_main_data_table$query_cluster,
                                                    as.character(actino789_main_data_table$ref_cluster),"different")

actino789_main_data_table$subcluster_compare <- ifelse(actino789_main_data_table$ref_subcluster==actino789_main_data_table$query_subcluster,
                                                       as.character(actino789_main_data_table$ref_subcluster),"different")

#Convert to factor class
actino789_main_data_table$host_genus_compare <- as.factor(actino789_main_data_table$host_genus_compare)
actino789_main_data_table$cluster_compare <- as.factor(actino789_main_data_table$cluster_compare)
actino789_main_data_table$subcluster_compare <- as.factor(actino789_main_data_table$subcluster_compare)

#Convert all Unspecified fields to NA missing value
actino789_main_data_table[actino789_main_data_table == "Unspecified"] <- NA



#Fig. 8a
plot_gcd_values_by_phage(actino789_main_data_table,'monty')

#Fig. 8a
plot_gcd_values_by_phage(actino789_main_data_table,'gma2')





###Gene content dissimilarity analysis for computing how isolated (discrete) each phage is

#Output data for python script
#Create an empty column to matching input data file formatting for downstream python script
actino789_main_data_table$empty <- 0
cyano209_main_data_table$empty <- 0

actino_output_subset <- subset(actino789_main_data_table,select=c('reference',
                                                              'query',
                                                              'empty',
                                                              'pham_pham_dissimilarity'))

cyano_output_subset <- subset(cyano209_main_data_table,select=c('reference',
                                                              'query',
                                                              'empty',
                                                              'pham_pham_dissimilarity'))

actino_cyano_output_subset <- rbind(actino_output_subset,cyano_output_subset)
write.table(actino_cyano_output_subset,"actino789_cyano209_pham_data.csv",sep=",",row.names = FALSE,col.names=FALSE,quote=FALSE)



#Import data from python script
#All gcd_gap data
gcd_gap_data <- read.csv("gcd_gap_analysis_all.csv",sep=",",header=TRUE)


#Fig. 8b
gcd_gap_data_sorted <- gcd_gap_data[order(gcd_gap_data$gcd_gap_max,decreasing=TRUE),]
par(mar=c(4,8,4,4))
plot(gcd_gap_data_sorted$gcd_gap_max,ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




###Analyze lineages (remove singletons)
gcd_gap_by_lineage_no_singleton <- subset(gcd_gap_data,gcd_gap_data$lineage != "Unspecified" &
                                          gcd_gap_data$lineage != "Singleton")

gcd_gap_by_lineage_no_singleton$lineage <- factor(gcd_gap_by_lineage_no_singleton$lineage)


#Fig. 8g
par(mar=c(8,12,2,4))
boxplot(gcd_gap_by_lineage_no_singleton$gcd_gap_max ~ gcd_gap_by_lineage_no_singleton$lineage,ylim=c(0,1),las=2,cex.axis=1,whisklty=0,staplelty=0,
        range=0,varwidth=FALSE,at=rank(tapply(gcd_gap_by_lineage_no_singleton$gcd_gap_max,gcd_gap_by_lineage_no_singleton$lineage,median)))
stripchart(gcd_gap_by_lineage_no_singleton$gcd_gap_max ~ gcd_gap_by_lineage_no_singleton$lineage,ylim=c(0,1),add=TRUE,vertical=TRUE,method="jitter",
           pch=19,cex=0.5,at=rank(tapply(gcd_gap_by_lineage_no_singleton$gcd_gap_max,gcd_gap_by_lineage_no_singleton$lineage,median)))







###Analyze Gordonia phages
gordonia_phages <- subset(gcd_gap_data,gcd_gap_data$host_genus == "Gordonia")
gordonia_phages$cluster <- factor(gordonia_phages$cluster)

#Fig. 8c
par(mar=c(8,12,2,4))
stripchart(gordonia_phages$gcd_gap_max ~ gordonia_phages$cluster,ylim=c(0,1),vertical=TRUE,method="jitter",pch=19,cex=1,
           las=2,cex.axis=1,at=rank(tapply(gordonia_phages$gcd_gap_max,gordonia_phages$cluster,median)))










#Analyze MaxGCDGap by Subcluster
#This data subset does not contain:
#1. Low frequency Subclusters
#2. Cyanobacteriophages
#3. Actinobacteriophages that have not been subclustered (such as Singletons)
gcd_gap_data_subcluster <- read.csv("gcd_gap_analysis_subcluster.csv",sep=",",header=TRUE)

#Fig. 8d
par(mar=c(8,12,2,4))
boxplot(gcd_gap_data_subcluster$gcd_gap_max ~ gcd_gap_data_subcluster$subcluster,ylim=c(0,1),las=2,cex.axis=1,whisklty=0,staplelty=0,range=0,
        varwidth=FALSE,at=rank(tapply(gcd_gap_data_subcluster$gcd_gap_max,gcd_gap_data_subcluster$subcluster,median)))
stripchart(gcd_gap_data_subcluster$gcd_gap_max ~ gcd_gap_data_subcluster$subcluster,ylim=c(0,1),add=TRUE,vertical=TRUE,method="jitter",
           pch=19,cex=0.5,at=rank(tapply(gcd_gap_data_subcluster$gcd_gap_max,gcd_gap_data_subcluster$subcluster,median)))



#Analyze MaxGCDGap by Cluster
#This data subset does not contain:
#1. Low frequency Clusters (although it does retain Singletons)
#2. Cyanobacteriophages
gcd_gap_data_cluster <- read.csv("gcd_gap_analysis_cluster.csv",sep=",",header=TRUE)



#Fig. 8e
par(mar=c(8,12,2,4))
boxplot(gcd_gap_data_cluster$gcd_gap_max ~ gcd_gap_data_cluster$cluster,ylim=c(0,1),las=2,cex.axis=1,whisklty=0,staplelty=0,range=0,
        varwidth=FALSE,at=rank(tapply(gcd_gap_data_cluster$gcd_gap_max,gcd_gap_data_cluster$cluster,median)))
stripchart(gcd_gap_data_cluster$gcd_gap_max ~ gcd_gap_data_cluster$cluster,ylim=c(0,1),add=TRUE,vertical=TRUE,method="jitter",
           pch=19,cex=0.5,at=rank(tapply(gcd_gap_data_cluster$gcd_gap_max,gcd_gap_data_cluster$cluster,median)))



#Analyze MaxGCDGap by host genus
#This data subset does not contain:
#1. Low frequency host genera
gcd_gap_data_genus <- read.csv("gcd_gap_analysis_genus.csv",sep=",",header=TRUE)

#Fig. 8f  
par(mar=c(8,12,2,4))
boxplot(gcd_gap_data_genus$gcd_gap_max ~ gcd_gap_data_genus$host_genus,ylim=c(0,1),las=2,cex.axis=1,whisklty=0,staplelty=0,range=0,
        varwidth=FALSE,at=rank(tapply(gcd_gap_data_genus$gcd_gap_max,gcd_gap_data_genus$host_genus,median)))
stripchart(gcd_gap_data_genus$gcd_gap_max ~ gcd_gap_data_genus$host_genus,ylim=c(0,1),add=TRUE,vertical=TRUE,method="jitter",
           pch=19,cex=0.5,at=rank(tapply(gcd_gap_data_genus$gcd_gap_max,gcd_gap_data_genus$host_genus,median)))














