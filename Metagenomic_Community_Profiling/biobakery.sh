#######Knead Data#######

###build database from Neotropical primate genomes available at NCBI Saimiri boliviensis boliviensis – GCA\_000235385.1; Callithrix jacchus – GCA\_002754865.1; Aotus nancymaae – GCA\_000952055.2; Cebus capucinus imitator - GCA\_001604975.1
### Knead Data v0.7.10
### Bowtie2 v 2.4.2 
### Trimmomatic-0.36


mkdir biobakery
cd    biobakery
mkdir biobakery/kneaddata_db #main directory for kneaddata analysis 
mkdir biobakery/kneaddata_db/kneaddata_output


bowtie2-build NeotropicalPrimate.fasta NeotropicalPrimate_db

#filter input data, input data has to be in the same folder as the database it seems 


kneaddata --input sample.R1.fastq.gz --input sample.R2.fastq.gz -db NeotropicalPrimate_db --output biobakery/kneaddata_db/kneaddata_output --trimmomatic  /user/bin/Trimmomatic-0.39/


#######Metaphlan#######
### Metaphlan v3.0.7
 
#create conda environment
conda create --name mpa -c bioconda python=3.7 metaphlan
conda activate mpa

#install metaphlan database
metaphlan --install --index mpa_v30_CHOCOPhlAn_201901  --bowtie2db metaphlan_db # location of database biobakery/metaphlan_db
 
 
#run metaphalan to estimate absolute read counts for bacterial taxon 
mkdir metaphlan_db #main directory for metaphalan directory
mkdir metaphlan_db/metaphlan_output

metaphlan kneaddata_db/kneaddata_output/sample.R1_kneaddata_paired_1.fastq,kneaddata_db/kneaddata_output/sample.R1_kneaddata_paired_2.fastq --bowtie2out metaphlan_db/metaphlan_output/sample.bowtie2.bz2 --input_type fastq -o metaphlan_db/metaphlan_output/sample.profiled_metagenome.txt -t  rel_ab_w_read_stat



#######Humann #######
### HUMAnN3  v3.0.0.alpha.4
#pair-ended files for each sample were concatenated together as input for humann

mkdir humann_db #main directory for humann analysis 
mkdir humann_db/input  #copy concatenated fastqs here 
mkdir humann_db/out #humann outputs go here 

cd humann_db
#download databases 
humann_databases --download chocophlan full humann_db

humann_databases --download uniref uniref90_diamond humann_db

#run humann
humann --input sample.concat.fastq.gz --output humann_db/out/sample_fastq

#make directory into which humann gene family outputs will be copied to merge them into a single table 
#mkdir humann_db/gene_fams
#cp humann_db/out/*/*genefamilies.tsv  humann_db/gene_fams/

#join humann output tables together   
#humann_join_tables -i gene_fams -o bacterioma_genefamilies.tsv --file_name genefamilies
  
# Normalizing RPKs to relative abundance
#mv bacterioma_genefamilies.tsv gene_fams/
#cd gene_fams

#humann_renorm_table -i bacterioma_genefamilies.tsv -o bacterioma_genefamilies-cpm.tsv --units cpm  --update-snames
   
   
# Regrouping genes to other functional categories
#humann_regroup_table --input bacterioma_genefamilies-cpm.tsv --output bacterioma_genefamilies.rxn-cpm.tsv --groups uniref90_rxn
 
 
#Attaching names to metacyc features
#humann_rename_table --input bacterioma_genefamilies.rxn-cpm.tsv  --output  bacterioma_genefamilies.metacyc-cpm.tsv  --names metacyc-rxn

 


### pathway abundances 
mkdir pathways 
mv *pathabundance.tsv pathways



for f in *out/*pathabundance.tsv
do 
cp $f pathways
done 


cd pathways
 
# Normalizing RPKs to relative abundance
for f in *tsv
do
echo $f 
humann_renorm_table -i $f -o ${f%tsv}cmp.tsv --units cpm  --update-snames
done 


mkdir normalized 
mv *cmp.tsv normalized

#join normalized table 
humann_join_tables -i normalized/ -o bacterioma_pathways.tsv

mv bacterioma_pathways.tsv normalized




#plot out significant pathways from humann analysis based on  results from  MaAsLin 2.022
for f in PWY-5154 FERMENTATION-PWY GLYCOLYSIS-E-D HEME-BIOSYNTHESIS-II HEXITOLDEGSUPER-PWY P161-PWY P4-PWY P461-PWY PHOSLIPSYN-PWY PWY-5791 PWY-5837 PWY-5840 PWY-5897 PWY-5898 PWY-5899 PWY-5659 ARG+POLYAMINE-SYN COA-PWY-1 FAO-PWY HISDEG-PWY POLYAMSYN-PWY PWY-1269 PWY-5083 PWY-5136 PWY-5173 PWY-5188 PWY-5345 PWY-5675 PWY-5918 PWY-6147_6 PWY-6163 PWY-621 PWY-4242 NAGLIPASYN-PWY ARO-PWY COA-PWY PWY-5138 PWY-5100
do
echo $f
humann_barplot --input bacterioma_pathways_sample.info.tsv --output ${f}.png --focal-feature $f --sort sum metadata --focal-metadatum  ENV
done


#values of pathways that didnt work in for loop 
humann_barplot --input bacterioma_pathways_sample.info.tsv --output PWY-6147.png --focal-feature PWY-6147 --sort sum metadata --focal-metadatum  ENV

for f in PWY-5154 FERMENTATION-PWY GLYCOLYSIS-E-D HEME-BIOSYNTHESIS-II HEXITOLDEGSUPER-PWY P161-PWY P4-PWY P461-PWY PHOSLIPSYN-PWY PWY-5791 PWY-5837 PWY-5840 PWY-5897 PWY-5898 PWY-5899 PWY-5659 ARG+POLYAMINE-SYN COA-PWY-1 FAO-PWY HISDEG-PWY POLYAMSYN-PWY PWY-1269 PWY-5083 PWY-5136 PWY-5173 PWY-5188 PWY-5345 PWY-5675 PWY-5918 PWY-6147_6 PWY-6163 PWY-621 PWY-4242 NAGLIPASYN-PWY ARO-PWY COA-PWY PWY-5138 PWY-5100; do humann_barplot --input bacterioma_pathways_sample.info.tsv --output ${f}.png --focal-feature $f --sort sum metadata --focal-metadatum  ENV -a normalize --remove-zeroes; done


#values of pathways that didnt work in for loop 
humann_barplot --input bacterioma_pathways_sample.info.tsv --output PWY-6147.png --focal-feature $f --sort sum metadata --focal-metadatum  ENV -a normalize --remove-zeroes



