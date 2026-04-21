require(dplyr)
require(ggplot2)
require(phyloseq)
require(microViz)
require(CompStellation)
require(vegan)
#read in metadata table
metadata_antony<-read.delim("C:/Users/kmmuf/Downloads/metadata 1.tsv")
#samples are rows, variables are columns
#using enterosignatures
#samples are rows, variables are columns
enter_antony<-read.delim("C:/Users/kmmuf/Downloads/aligned_model_samples.tsv")
#vv run this if samples are columns vv
#tenter_antony<-t(enter_antony)

#force to become a df
tenter_antony<-as.data.frame(enter_antony)
#making sure that there is a matching "internal_id" to join across
tenter_antony$internal_id<-tenter_antony$X
list(tenter_antony$internal_id %in% metadata_antony$internal_id)
#joining to produce one metadata+data table
anthonys_df<-left_join(metadata_antony,tenter_antony,join_by("internal_id"))
#removing some extraneous columns
anthonys_df_minimal<-anthonys_df[,c(2,17, 19,41:47)]
anthonys_df_minimalG<-subset(anthonys_df_minimal, !is.na(anthonys_df_minimal$gender))
colnames(anthonys_df_minimalG)<-c("internal_id","gender","geographic_location", "diet","X","S1","S2","S3","S4","S5")
#subset to variables of interest
anthonys_df_minimalG_sub<-subset(anthonys_df_minimalG,geographic_location %in% c("Fiji", "USA","China","Cameroon", "Madagascar", "Ethiopia", "Ireland", "United Kingdom", "Indonesia","Singapore","Tanzania", "Peru")  )
#this is an example of a non-phyloseq object being used for plot_taxa_star
a<-plot_taxa_star(
  physeq = anthonys_df_minimalG_sub,
  taxa_names = c("S1","S2","S3","S4","S5"),
  sample_var = "geographic_location",
  samplecolumn = "internal_id",
  error_bar = "SE",
  taxa_rank = "OTU",
  plot_order = "hclust")+ggtitle("Enterosignatures across countries")
ggsave("C:/Users/kmmuf/Downloads/enterotypes.png", a, width = 8, height = 6, dpi = 300, units = "in")

####### Build physeq obj for melissa ########
#load otu table
anthony_genera<-read.delim("C:/Users/kmmuf/Downloads/x.tsv")
#designate sample dataframe as sample frame
metadata_antony_sam<-sample_data(metadata_antony)
metadata_antony_sam@row.names<-metadata_antony_sam$internal_id

#making sure each otu has a unique id
rownames(anthony_genera)<-anthony_genera[,1]
#removing the spurious original genus name from the otu table
anthony_genera<-anthony_genera[,-1]
#make teh row names into a separate taxa matrix
taxa_TABLE<-as.matrix(rownames(anthony_genera))
#make sure otus in "taxa" and "otu" have same name
rownames(taxa_TABLE)<-rownames(anthony_genera)
#force it to be a taxa table
taxa_tables<-tax_table(taxa_TABLE)
#make ps object
physeq <- phyloseq(
  otu_table(anthony_genera, taxa_are_rows = TRUE),tax_table(taxa_tables),
  sample_data(metadata_antony_sam))

#subset to samples of interest
antgen<-subset_samples(physeq, geographic_location %in% c("Fiji", "USA","China","Cameroon", "Madagascar", "Ethiopia", "Ireland", "United Kingdom", "Indonesia","Singapore","Tanzania", "Peru") )
#remove extraneous zeros
antgen<-subset_taxa(antgen, taxa_sums(antgen)> 0)
#ordinate
ord<-ordinate(antgen, method = "PCoA",distance = "bray")
#plot beta diversity
plot_ord_comar<-plot_ordination(antgen,ord,color = "geographic_location")+geom_point(size = 2.5)
#plot alpha
plot_richness(antgen)

#extra nugget, use kaden's pcoa plot
anything<-plot_ordi_star(
  physeq = antgen,
  sample_var = "geographic_location",
  distance = "bray",
  view_type = "separate",
  plot_order = "hclust",
  fill_alpha = 0.3
)
anything2<-anything+ ggtitle("PCoA of genus level gut communities by country")+theme(legend.position = "none")
ggsave("C:/Users/kmmuf/Downloads/genera_anthony.png", anything2, width = 11, height = 11, dpi = 300, units = "in")
