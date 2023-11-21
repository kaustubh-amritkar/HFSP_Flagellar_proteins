## This file will have the bash commanda used during this project work

# The command for obtaining genome accession ids from protein accession ids
## Location (The output is written in): `/mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/`
for i in `cat ../data/sequences/protein_list.txt`; do printf $i","; esearch -db protein -query $i | elink -target nuccore | elink -target assembly | esummary | xtract -pattern AssemblyAccession -element AssemblyAccession; done > ../data/sequences/prot_genome_accession_mapping.csv

# Command for downloading the genome prot sequence from the ftp link
## Location: `/mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/pH_prediction/genome_datset`
while read -r line; do var=`echo $line | cut -d"," -f3`; add=`echo $var | rev | cut -d"/" -f1 | rev`; wget $var"/"$add"_protein.faa.gz"; done < ../../Caroline_data/species_list_30082023_HFSP.csv
# Unzip the .faa.gz genome fasta
for i in `ls`; do gzip -d $i; done

#Created file pfam_list_pH_genes.txt based on the list of genes of interested in the .ipynb code

# pH_hmms is the file with hmm profile for each of the genes of interest, provided by Morgon
# Command for performing the hmmer search for the pfam profiles against the genomes
## Location: `/mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/pH_prediction`
for file in `ls ./genome_datset/*.faa`; do hmmsearch --tblout $file.tbl.out -T 10 --domT 10 --incT 10 --incdomT 10 --noali pH_hmms $file; done

for i in `ls *.tbl.out`; do var=`echo $i | rev | cut -d"." -f2- | rev`; sed "/^#/d" $i | cut -c33-53 | sort -n | uniq > $var".pFam.pres.txt"; done

# The following command obtains the presence/absence info for each of the pfam onto the genome (open the output file and clean the first row)
for file in `ls ./hmmsearch_on_genome_results/*.faa.tbl.pFam.pres.txt`; do for pfam in `cat pfam_list_pH_genes.txt`; do genome_id=`echo $file | sed "s/_protein.faa.tbl.pFam.pres.txt//g" | cut -d"/" -f3`; value=`grep $pfam $file | wc -l`; echo $genome_id","$pfam","$value; done; done > ind_pfam_gene_presence_absence_on_genome.csv


for i in `cat ../sequences/prot_genome_accession_mapping.csv`; do phy_prot=`echo $i | cut -d"," -f1`; phy_genome=`echo $i | cut -d"," -f2`; phy_genome_num=`echo $phy_genome | sed "s/GCF_//g"`; genome_id=`grep $phy_genome_num ./Caroline_genomes_pH_predictions.csv | cut -d"," -f2`; ph_val=`grep $phy_genome_num ./Caroline_genomes_pH_predictions.csv | rev | cut -d"," -f1 | rev`; species_info=`grep $genome_id ../Caroline_data/species_list_30082023_HFSP.csv | cut -d"," -f4-`; echo $phy_prot","$phy_genome","$ph_val","$species_info; done > phylogeny_features_info.csv

# Command for obtaining the organism name for the protein acc identifiers
for i in `cat ./prot_genome_accession_mapping.csv`; do prot_acc=`echo $i | cut -d"," -f1`; organism=`esearch -db protein -query $prot_acc | elink -target nuccore | elink -target assembly | esummary | xtract -pattern Organism -element Organism`; echo $i","$organism; done > prot_genome_accession_organism_mapping.csv
sed "s/(.*)//g" prot_genome_accession_organism_mapping.csv > prot_genome_accession_organism_name_mapping.csv
sed -i "s/ /_/g" prot_genome_accession_organism_name_mapping.csv
sed -i "s/_$//g" prot_genome_accession_organism_name_mapping.csv
sed -i "s/\//_/g" prot_genome_accession_organism_name_mapping.csv

cut -d"," -f2,3 prot_genome_accession_organism_name_mapping.csv | sort -n | uniq | tail -n +2 > uniq_genome_accession_organism_name_mapping.csv

### OGT prediction work from here ###
# Downloading the nucleotide genomes for the genomic identifiers for the phylogeny protein entries
for i in `cut -d"," -f2 prot_genome_accession_organism_name_mapping.csv | sort -n | uniq`; do esearch -db assembly -query $i | elink -target nucleotide -name \ assembly_nuccore_refseq | efetch -format fasta > ./nuc_genome/$i.fa; done

# For creating the genomes directory for OGT prediction
## Location: `/mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/analysis/OGT_prediction/prediction_HFSP`
for i in `cut -d"," -f2,3 /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/prot_genome_accession_organism_name_mapping.csv | sort -n | uniq | tail -n +2`; do genome_acc=`echo $i | cut -d"," -f1`; org_name=`echo $i | cut -d"," -f2`; echo $org_name; mkdir genomes/$org_name; cp /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/nuc_genome/$genome_acc.fa ./genomes/$org_name/; done
##for i in `cut -d"," -f2,3 /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/prot_genome_accession_organism_name_mapping.csv | sort -n | uniq | tail -n +2`; do genome_acc=`echo $i | cut -d"," -f1`; genome_gca_id=`echo $genome_acc | sed "s/GCF/GCA/g"`; org_name=`echo $i | cut -d"," -f2`; echo $org_name; mkdir genomes/$org_name; cp /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/pH_prediction/genome_datset/$genome_gca_id*.faa ./genomes/$org_name/; done
find ./ -type d -empty ##For finding the directories that are empty; deleted the directories corresponding to these
cut -d"," -f2,3 prot_genome_accession_organism_name_mapping.csv | sort -n | uniq > uniq_genome_accession_organism_name_mapping.csv  ##Made a different file with these unique entries; then manually removed the entries that don't have a corresponding genome seq

# Run the ogt_prediction_preparing_taxa_file.ipynb jupyter notebook to create the `species_taxanomic.txt` file; Manually edit the first header to `species`

# gzip the genomic fasta files, to run the OGT-predictor
## Location: `/mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/analysis/OGT_prediction/prediction_HFSP/genomes/`
for i in `ls */*.faa`; do gzip $i; done

# Making the genomes_retrived.txt file
for i in `cut -d"," -f2,3 /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/prot_genome_accession_organism_name_mapping.csv | sort -n | uniq | tail -n +2`; do genome_acc=`echo $i | cut -d"," -f1`; org_name=`echo $i | cut -d"," -f2`; genome_seqs_file=`ls /mnt/researchdrive/Kaustubh/HFSP-Flagellar_proteins/data/sequences/nuc_gen
ome/$genome_acc.fa | rev | cut -d"/" -f1 | rev`; echo -e $genome_seqs_file"\t"$org_name; done > genomes_retrieved.txt
sed -i "s/.fa/.fa.gz/g" genomes_retrieved.txt


# For obtaining the predicted ogt values based on the common entries from experimental dataset
for i in `cat genome_exp_ogt_data_common_with_HFSP.csv`; do genome_num=`echo $i | cut -d"," -f1 | cut -d"_" -f2`; var=`echo "GCF_"$genome_num`; org_name=`grep $var ../prediction_HFSP/genomes_retrieved_wo_sus.txt | cut -f2`; pred_ogt=`grep $org_name ../prediction_HFSP/newly_predicted_OGTs.txt | cut -f2`; echo $var","$org_name","$i","$pred_ogt; done > ogt_exp_pred_comparison.csv


# Combining the results from different sections
for i in `cat ../pH_prediction/phylogeny_labels_with_pH.csv | tail -n +2`; do phy_label=`echo $i | cut -d"," -f1`; prot_id=`echo $i | cut -d"," -f2`; gcf_id=`echo $i | cut -d"," -f3`; genome_num=`echo $i | cut -d"," -f3 | cut -d"_" -f2 | cut -d"." -f1`; pH_pred=`echo $i | cut -d"," -f4`; ogt_map=`grep $gcf_id ../ogt_prediction/genome_org_taxa_ogt_pred.csv`; org_name=`echo $ogt_map | cut -d"," -f2`; lineage=`echo $ogt_map | cut -d"," -f3-7`; ogt_pred=`echo $ogt_map | cut -d"," -f8`; species_info=`grep $gcf_id ../Caroline_data/species_list_30082023_HFSP.csv`; flagellum_info=`echo $species_info | cut -d"," -f5`; motility_info=`echo $species_info | cut -d"," -f6`; echo $phy_label","$prot_id","$gcf_id","$lineage","$pH_pred","$ogt_pred","$flagellum_info","$motility_info; done > phy_genome_lineage_pH_ogt_species-info.csv
## Some manual editing after that


###############
### THE FOLLOWING IS FOR THE WORK RELATED TO SPECIES TREE & MAPPING THE GENOMIC FEATURES ONTO THE SPECIES TREE ###
###############

# Downloading the protein sequence genomes
while read -r line; do assembly_id=`echo $line | cut -d"," -f1`; GCA_id=`echo $assembly_id | cut -d"_" -f2 | cut -d"." -f1`; x1=`echo $GCA_id | cut -c1-3`; x2=`echo $GCA_id | cut -c4-6`; x3=`echo $GCA_id | cut -c7-9`; echo $assembly_id; wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"$x1"/"$x2"/"$x3"/"$assembly_id"/"$assembly_id"_protein.faa.gz"; done < ../bactgenomes.csv

# Downloading the nucleotide sequence genomes
while read -r line; do assembly_id=`echo $line | cut -d"," -f1`; GCA_id=`echo $assembly_id | cut -d"_" -f2 | cut -d"." -f1`; x1=`echo $GCA_id | cut -c1-3`; x2=`echo $GCA_id | cut -c4-6`; x3=`echo $GCA_id | cut -c7-9`; echo $assembly_id; wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"$x1"/"$x2"/"$x3"/"$assembly_id"/"$assembly_id"_genomic.fna.gz"; done < ../bactgenomes.csv

# Shaping the files for running the OGT-prediction code like the genomes_retrieved.txt & the genomes folder
cut -d"," -f1,7 bactgenomes.csv > genomes_retrieved.txt     #Remove the header manually
sed -i "s/\//_/g" genomes_retrieved.txt
sed -i "s/_(.*//g" genomes_retrieved.txt
sed -i "s/,/_genomic.fna.gz\t/g" genomes_retrieved.txt


# Command for replacing the horrible phy_label with the better assembly_id as the phy_label
### `16srna_species_tree.tree` file was made by saving the figtree output
for i in `cut -f1 genomes_retrieved_req.txt | sed "s/_genomic.fna.gz//g"`; do phy_label=`grep $i species_tree_headers.txt`; phy_label_upd=$(sed 's/[\/&]/\\&/g' <<< $phy_label); sed -i "s/$phy_label_upd/$i/g" 16srna_species_tree.tree; done
sed -i "s/'//g" 16srna_species_tree.tree

