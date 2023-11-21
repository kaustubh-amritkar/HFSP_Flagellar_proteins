library(ggtree)
library(ape)
library(Biostrings)
library(ggimage)
library(ggnewscale)
library(svglite)
library(viridis)
library(phytools)
library(randomcoloR)

MotA = ape::read.tree(file = "/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/data/Caroline_data/motA444_tree_uncolored.treefile")

Mota_features_df = data.frame(read.csv("/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/data/combine_ph_ogt_taxa_festures/phy_genome_lineage_pH_ogt_species-info-2.csv"))
Mota_features_df[Mota_features_df == ""] <- NA

pH_pred = Mota_features_df[9]
rownames(pH_pred) = Mota_features_df[,1]

OGT_pred = Mota_features_df[10]
rownames(OGT_pred) = Mota_features_df[,1]

phylum = Mota_features_df[5]
rownames(phylum) = Mota_features_df[,1]

Flagellum = Mota_features_df[11]
rownames(Flagellum) = Mota_features_df[,1]

Motility = Mota_features_df[12]
rownames(Motility) = Mota_features_df[,1]

prot_type = Mota_features_df[13]
rownames(prot_type) = Mota_features_df[,1]

MotA_tree_branch = ggtree(MotA, size=0.2, aes()) %<+% Mota_features_df + geom_tree(aes(color=prot_type), size=1.0) + geom_treescale()
ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/MotA_tree_with_prot_type_branch.svg", plot = MotA_tree_branch, width = 20, height = 16)

MotA_tree = ggtree(MotA, size=0.2, aes()) %<+% Mota_features_df + geom_tippoint(aes(color=phylum), size=2.5) + geom_treescale()

MotA_tree_1 = gheatmap(MotA_tree, pH_pred, low="orange", high="#008080", offset=0.5, width=0.1, colnames=FALSE) + scale_x_ggtree() + scale_y_continuous(limits=c(4,8), expand = c(0, 0.3))
MotA_tree_2 = MotA_tree_1 + new_scale_fill()
MotA_tree_3 = gheatmap(MotA_tree_2, OGT_pred, low="blue", high="red", offset=1.0, width=0.1, colnames=FALSE) + scale_x_ggtree() + scale_y_continuous(limits=c(20,90), expand = c(0, 0.3))
MotA_tree_4 = MotA_tree_3 + new_scale_fill()
MotA_tree_5 = gheatmap(MotA_tree_4, Flagellum, offset = 1.5, width = 0.1, colnames = FALSE) + scale_fill_manual(breaks=c("Yes","Yes?","No","?","Pilus?"), values=c("steelblue","cyan","brown","yellow","tan"), na.value = NA)
MotA_tree_6 = MotA_tree_5 + new_scale_fill()
MotA_tree_7 = gheatmap(MotA_tree_6, Motility, offset = 2.0, width = 0.1, colnames = FALSE) + scale_fill_manual(breaks=c("Gliding","Motile","Motile?","Nonmotile"), values=c("steelblue","brown","yellow","tan"), na.value = NA)

# MotA_tree_5 = gheatmap(MotA_tree_4, Flagellum, offset = 1.5, width = 0.1, colnames = FALSE) + scale_fill_manual(breaks=c("Yes","Yes?","No","?","Pilus?"), values=distinctColorPalette(5), na.value = NA)
# MotA_tree_7 = gheatmap(MotA_tree_6, Motility, offset = 2.0, width = 0.1, colnames = FALSE) + scale_fill_manual(breaks=c("Gliding","Motile","Motile?","Nonmotile"), values=distinctColorPalette(4), na.value = NA)

ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/MotA_tree_with_taxa_env_features.svg", plot = MotA_tree_7, width = 20, height = 16)

Mota_pH_pred_distribution = ggplot(Mota_features_df, aes(x=pH_pred)) + geom_histogram(color="black", fill="darkgray") + theme_bw()
ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/MotA_pH_pred_distribution.png", plot = Mota_pH_pred_distribution)

Mota_OGT_pred_distribution = ggplot(Mota_features_df, aes(x=OGT_pred)) + geom_histogram(color="black", fill="darkgray") + theme_bw()
ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/MotA_OGT_pred_distribution.png", plot = Mota_OGT_pred_distribution)

Mota_pH_OGT_scatter_plot = ggplot(Mota_features_df, aes(x=OGT_pred, y=pH_pred, color=phylum)) + geom_point()
ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/MotA_pH_OGT_scatter_plot.png", plot = Mota_pH_OGT_scatter_plot)

### The Following Part is for the Species Tree mapping  ###
species = ape::read.tree(file = "/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/data/species_tree/16srna_species_tree.tree")

species_df = data.frame(read.csv("/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/data/species_tree/species_tree_genomic_predictions_for_mapping_short_phy_label.csv"))
species_df[species_df == ""] <- NA

sp_pH_pred = species_df[3]
rownames(sp_pH_pred) = species_df[,1]

sp_ogt_pred = species_df[4]
rownames(sp_ogt_pred) = species_df[,1]

species_tree = ggtree(species, size=0.2, aes()) %<+% species_df + geom_tippoint(aes(color=phylum), size=1.0) + geom_treescale()

species_tree_1 = gheatmap(species_tree, sp_ogt_pred, low="blue", high="red", offset=0.05, width=0.1, colnames=FALSE)

species_tree_2 = species_tree_1 + new_scale_fill()
species_tree_3 = gheatmap(species_tree_2, sp_pH_pred, low="orange", high="#008080", offset=0.15, width=0.1, colnames=FALSE)

ggsave(file="/Volumes/bkacar/Kaustubh/HFSP-Flagellar_proteins/figures/Species_tree_with_taxa_env_features.svg", plot = species_tree_3, width = 20, height = 16)

ggplot(species_df, aes(x=sp_pH_pred)) + geom_histogram(color="black", fill="darkgray") + theme_bw()

ggplot(species_df, aes(x=sp_ogt_pred)) + geom_histogram(color="black", fill="darkgray") + theme_bw()
