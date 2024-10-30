#!/bin/bash

#####################################################################
#
#########     GBS_filter_GS_v19.2.sh    17 sept 2024    #############
#
#####################################################################

#SBATCH -J GBfi-v19
#SBATCH -o %j.log
#SBATCH -c 4
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vtbos@ulaval.ca
#SBATCH --time=3:00:00
#SBATCH --mem=15G

#######  Variables to define ######
input_vcf_file=GAUP2023-BRTS2023_raw_renamed.vcf
path=/home/vtbos/GAPP/1-GAUP2023_BRTS2023/Prj-4/fast-gbs_v3/results	#path for your input_vcf_file
outputs_folder=/home/vtbos/GAPP/1-GAUP2023_BRTS2023/Prj-4/fast-gbs_v3/results/redo1-filtrations-includingROAP
scaffold_name=scaffold				# enter the name used in your .vcf file for scaffolds name in the CHROM column (ex: SCAFF, scaffold, ChrUn, scaf)
missingness_threshold=0.8			# enter the threshold value for missing data per individual
empty_word=EMPTY				# enter the character patern given to your empty wells/samples in the plates (ex: empty, EMPTY, NA, vide)
###################################

# loading the softwares
module load vcftools/0.1.16
module load tassel/5.2.29
module load java/jdk/1.8.0_102

# creating the directory for all outputs
mkdir $outputs_folder
cd $outputs_folder
name=${input_vcf_file::-4}

# defining the function for Tassel summaries
function tassel_summaries {
	run_pipeline.pl -Xms8g -Xmx10g -fork1 -vcf $name".vcf" -genotypeSummary site -export $PREFIX$EXTENTION"_Tassel_site_report" -runfork1
	run_pipeline.pl -Xms8g -Xmx10g -fork1 -vcf $name".vcf" -genotypeSummary taxa -export $PREFIX$EXTENTION"_Tassel_taxa_report" -runfork1
	run_pipeline.pl -Xms8g -Xmx10g -fork1 -vcf $name".vcf" -genotypeSummary overall -export $PREFIX$EXTENTION"_Tassel_oall_report" -runfork1
}

#########################################################
# 	FILTRATION
#########################################################

#removing the scaffolds from the VCF file
EXTENTION=_noScaf

grep -v $scaffold_name $path"/"$name".vcf" > ./$name$EXTENTION".vcf"
name=$name$EXTENTION

#1
EXTENTION=_onlyPASS

vcftools --vcf $name".vcf" --remove-filtered-all --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"

#2
EXTENTION=_rmIndels

vcftools --vcf $name".vcf" --remove-indels --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"
PREFIX=1
tassel_summaries

#
EXTENTION=_NoEmpty

grep "#CHROM" $path"/"$input_vcf_file | sed 's/\t/\n/g' | grep $empty_word > $outputs_folder"/empty_wells_list.txt"
vcftools --vcf $name".vcf" --remove $outputs_folder"/empty_wells_list.txt" --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"

#3
EXTENTION=_minMaxAll2

vcftools --vcf $name".vcf" --min-alleles 2 --max-alleles 2 --recode --out $name$EXTENTION
name=$name$EXTENTION
name2=$name
mv $name".recode.vcf" $name".vcf"
PREFIX=2
tassel_summaries


################################################################################
#	 LOOP 1
################################################################################

#beginning the loop for exclusion of lines with too many misssing data
EXTENTION=_mac4_mm02
vcftools --vcf $name".vcf" --mac 4 --max-missing 0.2 --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"

vcftools --vcf $name".vcf" --missing-indv
mv out.imiss $name".imiss"
name=$name".imiss"

#writting the variable "name' and 'threshold' so it can be inputed by the R script created below
echo $name > temp_input_name_for_R.txt
echo $missingness_threshold > threshold_for_missingness_per_individual.txt
echo $empty_word > character_string_for_empty_wells.txt

#The next lines creates a R script that outputs a text file with the name of the lines(individuals) that are over the threshold (this variable can be changed in the R script below) of missing data
echo '#!/usr/bin/env Rscript

# In this R script we will identify the outliers using the percentiles method for missingness per line on the VCF generated .imiss file

temp <- readLines("temp_input_name_for_R.txt")
data <- read.table(temp, header =TRUE)
empty_word <- readLines("character_string_for_empty_wells.txt")
threshold <- readLines("threshold_for_missingness_per_individual.txt")
threshold <- as.numeric(threshold)

data$content <- rep("DNA", nrow(data))
emptyLines <- grep(empty_word, data$INDV)
data[emptyLines, "content"] <- empty_word

#subset the data frame keeping only lines that content DNA (eclude empty wells)
data_noEmpty <- data[data$content == "DNA",]

#determine the outlier threshold
s <- summary(data_noEmpty$F_MISS)
q1 <- s[2]
q3 <- s[5]

IQR3 <- q3 + 3*(q3-q1)
IQR1.5 <- q3 + 1.5*(q3-q1)

outlier_threshold <- threshold

# tag the lines over the threshold
data$over_threshold <- rep("NO", nrow(data))
data[data$F_MISS >= outlier_threshold, "over_threshold"] <- "YES"

# exporting
to_remove <- data[data$over_threshold == "YES", "INDV"]
write.table(x = to_remove, file = "lines_removed_over_threshold_of_missing_data.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#plotting
library(ggplot2)
ggplot(data = data, mapping = aes(y = F_MISS, x ="")) +
  geom_jitter(mapping = aes(colour = content, shape = over_threshold), height = 0) +
  geom_hline(yintercept = c(outlier_threshold, IQR3, IQR1.5), show.legend = TRUE) +
  ggtitle("jitter_outliers_missing_genotypes") +
  xlab("Count") + ylab("Proportion missing per sample") +
  geom_text(aes(0,outlier_threshold, label = paste(outlier_threshold, "(manual threshold)"), vjust = -1, hjust = -0.1)) +
  geom_text(aes(0.5,IQR3, label = paste(IQR3,"(Q3 + 3IQR)"), vjust = -1, hjust = -0.1)) +
  geom_text(aes(0.5,IQR1.5, label =paste(IQR1.5, "(Q3 + 1.5IQR)"), vjust = -1, hjust = -0.1))
 
ggsave("jitter_outliers_missing_genotypes.pdf")

library(ggplot2)
ggplot(data, aes(x = F_MISS)) +
  geom_histogram() +
  theme_classic() +
  geom_vline(xintercept = c(outlier_threshold, IQR3, IQR1.5))+
  xlab("Proportion missing per sample") + ylab("Count") +
  geom_text(aes(outlier_threshold, outlier_threshold, label = paste(outlier_threshold, "(manual threshold)"), vjust = 0, hjust = -2, angle = 90)) +
  geom_text(aes(as.double(IQR3), as.double(IQR3), label = paste(IQR3,"(Q3 + 3IQR)"), vjust = 0, hjust = -2, angle = 90)) +
  geom_text(aes(as.double(IQR1.5), as.double(IQR1.5), label = paste(IQR1.5, "(Q3 + 1.5IQR)"), vjust = 0, hjust = -2, angle = 90)) +
  ggtitle("density_curve_outliers_missing_genotypes")
  
ggsave("density_curve_outliers_missing_genotypes.pdf")


' > Rscript_imiss_stats_to_remove_lines_over_missing_data_threshold.R

# then we run the R script    
Rscript Rscript_imiss_stats_to_remove_lines_over_missing_data_threshold.R

#rm temp_input_name_for_R.txt
#rm Rplots.pdf

# then we run filters in VCFtools
name=$name2
EXTENTION=_highNrm
vcftools --vcf $name".vcf" --remove lines_removed_over_threshold_of_missing_data.txt --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"
PREFIX=3
tassel_summaries
# this last line finished the loop to remove the lines with too much missing data 


EXTENTION=_mac4 
vcftools --vcf $name".vcf" --mac 4 --recode --out $name$EXTENTION
name=$name$EXTENTION
mv $name".recode.vcf" $name".vcf"
PREFIX=4
tassel_summaries

EXTENTION=_mm20
vcftools --vcf $name".vcf" --max-missing 0.2 --recode --out $name$EXTENTION
name=$name$EXTENTION
name3=$name
mv $name".recode.vcf" $name".vcf"
PREFIX=5
tassel_summaries

###########################################################################################
#	LOOP 2
###########################################################################################

# beginning of the loop to excludes sites with too high heterozygosity rate
# imputation with Beagle
EXTENTION=_imputed
java -Xmx25000m -jar /prg/beagle/5.0/beagle.jar gt=$name".vcf" out=$name$EXTENTION
name=$name$EXTENTION
gunzip $name".vcf.gz"

rm $name".log"

# creating the GT table to perform stats on heterozygousity per sites in R
vcftools --vcf $name'.vcf' --extract-FORMAT-info GT --out $name

sed 's/|/\//g' $name".GT.FORMAT" > $name".GT.FORMAT.t"

# writing the name of the file that will be the input for the R script
echo $name".GT.FORMAT.t" > hetero_site_filt_input_name.txt


# writing the R script
echo '#!/usr/bin/env Rscript

# R script to find outliers for heterozygosity per site
# produces an output file of sites to be removed that can be the input for vcftools (vcftools function; --exclude-positions)

temp2 <- readLines("hetero_site_filt_input_name.txt")
data <- read.table(temp2, header = T)

matrix <- as.matrix(data[ ,3:ncol(data)])
matrix <- gsub("0/1", "1/0", matrix)
matrix <- matrix == "1/0"
matrix <- matrix*1
Proportion.Heterozygous <- rep(NA, nrow(matrix))

for (i in 1:nrow(matrix)) {
  Proportion.Heterozygous[i] <- sum(matrix[i,])/ncol(matrix)
}

data <- cbind(data, Proportion.Heterozygous)

#calculating the outliers based on the interquartile range
QT <- summary(data$Proportion.Heterozygous)[c(2, 5)]
QT <- as.vector(QT)
Q1 <- QT[1]
Q3 <- QT[2]
IQR3 <- Q3+(3*(Q3-Q1))
threshold <- Q3+(1.5*(Q3-Q1))

# filtering sites over Q3+1.5IQR
data_export <- data[data$Proportion.Heterozygous > threshold, c(1, 2)]
data <- data[, c("CHROM", "POS", "Proportion.Heterozygous")]

# producing the output file
write.table(data_export, file = "Heterozygosity_Outliers_sites_Q3plus1.5IQR.txt", col.names = F, row.names = F, quote = F)

# plotting
library(ggplot2)
total_n_sites <- length(data$Proportion.Heterozygous)
over1.5IQR <- data$Proportion.Heterozygous >= threshold
o <- which(over1.5IQR)
n_sites_over_1.5IQR <- length(o)
over3IQR <- data$Proportion.Heterozygous >= IQR3
p <- which(over3IQR)
n_sites_over_3IQR <- length(p) 
n_section1 <- total_n_sites-n_sites_over_1.5IQR
perct_section1 <- (n_section1*100)/total_n_sites
p1 <- round(perct_section1, digits = 1)
n_section3 <- n_sites_over_3IQR
perct_section3 <- (n_section3*100)/total_n_sites
p3 <- round(perct_section3, digits = 1)
n_section2 <- total_n_sites-n_section1-n_section3
perct_section2 <- (n_section2*100)/total_n_sites
p2 <- round(perct_section2, digits = 1)
threshold <- as.numeric(format(round(threshold, 2), nsmall = 2))
IQR3 <- as.numeric(format(round(IQR3, 2), nsmall = 2))
nb_snps_tot <- length(data$Proportion.Heterozygous)
nb_snps_lost_1.5IQR <- length(which(data$Proportion.Heterozygous > threshold))
nb_snps_lost_3IQR <- length(which(data$Proportion.Heterozygous > IQR3))
nb_snps_lost_0.1 <- length(which(data$Proportion.Heterozygous > 0.1))

# plotting
ggplot(data, aes(x = Proportion.Heterozygous)) +
  geom_density() +
  geom_vline(xintercept = threshold) +
  geom_vline(xintercept = IQR3) +
  geom_vline(xintercept = 0.1) +
  ggtitle("Distribution frequency of heterozygosity per SNP") +
  xlab("Proportion of heterozygosity") + ylab("Frequency") +
  geom_text(aes(threshold , 0, label = paste(threshold, "(Q3 + 1.5IQR)", ", filtering out", nb_snps_lost_1.5IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 90)) +
  geom_text(aes(IQR3, IQR3, label = paste(IQR3,"(Q3 + 3IQR)", ", filtering out", nb_snps_lost_3IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 90)) +
  geom_text(aes(0.1, 0.1, label = paste("0.1", ", filtering out", nb_snps_lost_0.1, "SNPs"), vjust = 0, hjust = -0.3, angle = 90)) +
  geom_text(aes(0.1, 0.1, label = paste("Total nb of sites =", nb_snps_tot), vjust = 0, hjust = -2, angle = 0))
  
ggsave("density_curve_proportion_heterozygosity_per_sites.pdf")

ggplot(data = data, mapping = aes(y = Proportion.Heterozygous, x ="")) +
  geom_jitter(shape = 1, width = 1, color = "grey") +
  geom_violin(color = "black", alpha = 0.2) +
  geom_hline(yintercept = c(threshold, IQR3, 0.1), show.legend = TRUE, color = "black") +
  ggtitle("Proportion heterozygosity per SNP") +
  xlab("Sites (SNPs)") + ylab("Proportion heterozygosity") +
  geom_text(aes(0.1 , threshold, label = paste(threshold, "(Q3 + 1.5IQR)", ", filtering out", nb_snps_lost_1.5IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 0), color = "black") +
  geom_text(aes(0.1, IQR3, label = paste(IQR3,"(Q3 + 3IQR)", ", filtering out", nb_snps_lost_3IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 0), color = "black") +
  geom_text(aes(0.1, 0.1, label = paste("0.1", ", filtering out", nb_snps_lost_0.1, "SNPs"), vjust = 0, hjust = -0.3, angle = 0), color = "black") +
  geom_text(aes(0.1, 0.9, label = paste("Total nb of sites =", nb_snps_tot), vjust = 0, hjust = -2, angle = 0), color = "black")

ggsave("jitter_proportion_heterozygous_per_sites.pdf")
' > Rscript_hetero_site_filtering.R


# running the R script
Rscript Rscript_hetero_site_filtering.R


# removing the site with vcftools
EXTENTION=_highHrm
vcftools --vcf $name3".vcf" --exclude-positions Heterozygosity_Outliers_sites_Q3plus1.5IQR.txt --recode --out $name3$EXTENTION
name=$name3$EXTENTION
mv $name".recode.vcf" $name".vcf"
PREFIX=6
tassel_summaries


# imputation with Beagle
EXTENTION=_imputed
java -Xmx25000m -jar /prg/beagle/5.0/beagle.jar gt=$name".vcf" out=$name$EXTENTION
name=$name$EXTENTION
gunzip $name".vcf.gz"
PREFIX=7
tassel_summaries

rm $name".log"
rm Rplots.pdf

###############################################################
#		SUMMARY
###############################################################

# script to create summary graphics for all filtration steps
# writing the R script
echo '#!/usr/bin/env Rscript

library(ggplot2)

overall_reports <- list.files(pattern = "Tassel_oall_report1")

# importing all Tassel overall reports to a unique data frame
overall_report_allvcfs <- NULL
for (i in overall_reports) {
	
  dat <- read.table(i, sep = "\t" ,header = T, fill = F, na.strings="", stringsAsFactors = FALSE, dec = ",")
  name <- substr(i, 3, nchar(i)-24)
  colnames(dat)[2] <- name
  overall_report_allvcfs <- cbind(overall_report_allvcfs, dat[, name])
  colnames(overall_report_allvcfs)[ncol(overall_report_allvcfs)] <- name
  rownames(overall_report_allvcfs) <- dat$Stat.Type
}

# some data wrangling
overall_report_allvcfs <- as.data.frame(overall_report_allvcfs)
colnames(overall_report_allvcfs) <- gsub(" ", "_", colnames(overall_report_allvcfs))
rownames(overall_report_allvcfs) <- gsub(" ", "_", rownames(overall_report_allvcfs))
overall_report_allvcfs <- t(overall_report_allvcfs)
d <- data.frame(rownames(overall_report_allvcfs))
overall_report_allvcfs <- cbind(overall_report_allvcfs, d)
overall_report_allvcfs <- overall_report_allvcfs[, c(ncol(overall_report_allvcfs), 1:(ncol(overall_report_allvcfs)-1))]
colnames(overall_report_allvcfs)[1] <- "filter"
filter_names_odered <- substr(overall_reports, 3, nchar(overall_reports)-24)
overall_report_allvcfs$filter <- factor(overall_report_allvcfs$filter, levels = filter_names_odered)
write.table(overall_report_allvcfs, file = "tassel_overall_report_after_each_filtration_steps.txt", sep = "\t", quote = F, row.names = F)

# dataframe type 2 for missing and non missing sites
df_number_not_missing <- overall_report_allvcfs[, c("filter", "Number_Not_Missing")]
df_number_not_missing$data_type <- rep("number_not_missing", nrow(df_number_not_missing))
colnames(df_number_not_missing)[2] <- "SNPs_number"
df_Number_Missing <- overall_report_allvcfs[, c("filter", "Number_Missing")]
df_Number_Missing$data_type <- rep("number_missing", nrow(df_Number_Missing))
colnames(df_Number_Missing)[2] <- "SNPs_number"
df_miss_not_miss <- rbind(df_Number_Missing, df_number_not_missing)

# plots for missing and non_missing sites
ggplot(data = df_miss_not_miss, mapping = aes(x = as.numeric(filter), y = SNPs_number, color = data_type)) +
  geom_point() +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Number of Sites (SNPs) after different filtering steps", x = "Filters")

ggsave("SNP_quantity_after_each_filters.pdf")

# plot for proportion heterozygous sites
ggplot(data = overall_report_allvcfs, 
       mapping = aes(x = as.numeric(filter), y = Proportion_Heterozygous)) +
  geom_point() +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion heterozygous Sites (SNPs) after different filtering steps", x = "Filters")

ggsave("Proportion_heterozygous_sites_after_each_filters.pdf")

############ importing the taxa summary by Tassel
taxa_reports <- list.files(pattern = "Tassel_taxa_report")

# importing all Tassel taxa reports to a unique data frame
taxa_report_allvcfs <- NULL
for (i in taxa_reports) {

dat <- read.table(i, sep = "\t" ,header = T, fill = F, na.strings="", stringsAsFactors = FALSE, dec = ",")  
  name <- substr(i, 3, nchar(i)-24)
  dat$filters <- rep(name, nrow(dat))
  taxa_report_allvcfs <- rbind(taxa_report_allvcfs, dat)
}

# some data wrangling
taxa_report_allvcfs$filters <- as.factor(taxa_report_allvcfs$filters)
filter_names_odered <- substr(taxa_reports, 3, nchar(taxa_reports)-24)
taxa_report_allvcfs$filters <- factor(taxa_report_allvcfs$filters, levels = filter_names_odered)
taxa_report_allvcfs$Taxa <- gsub("^*", "s", taxa_report_allvcfs$Taxa)

## ggplot graphs
ggplot(data = taxa_report_allvcfs, mapping =aes(x = as.numeric(filters), y = Proportion.Missing, group = as.numeric(filters))) +
  geom_boxplot(outlier.shape = 3) + 
  geom_jitter(width = 0.1, alpha = 0) +
  geom_text(aes(label = Taxa), position = position_jitter(width=0.4,height=0.2), size = 1.5) +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion missing sites (SNPs) per taxa after different filtering steps", x = "Filters")

ggsave("Proportion_missing_sites_per_taxa_after_each_filters.pdf")

ggplot(data = taxa_report_allvcfs, mapping =aes(x = as.numeric(filters), y = Proportion.Heterozygous, group = as.numeric(filters))) +
  geom_boxplot(outlier.shape = 3) + 
  geom_jitter(width = 0.1, alpha = 0) +
  geom_text(aes(label = Taxa), position = position_jitter(width=0.4,height=0.2), size = 1.5) +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion heterozygous sites (SNPs) per taxa after different filtering steps", x = "Filters")

ggsave("Proportion_heterozygous_sites_per_taxa_after_each_filters.pdf")

# outputting the table
write.table(taxa_report_allvcfs[,-c(4, 8, 9)], file = "tassel_taxa_report_after_each_filtration_steps.txt", quote = FALSE, col.names = T, row.names =  F, sep = "\t")
' > Rscript_tassel_overall_and_taxa_reports_for_all_filtration_steps.R

# running the R script
Rscript Rscript_tassel_overall_and_taxa_reports_for_all_filtration_steps.R

# creating the Markdown script to print a report for all stats perfomed in this SLURM

echo '
---
title: "FastGBS_filters_v19.1.Rmd"
author: ""
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, echo = FALSE}
# In this R script we will identify the outliers using the percentiles method for missingness per line on the VCF generated .imiss file
library(ggplot2)
temp <- readLines("temp_input_name_for_R.txt")
data <- read.table(temp, header =TRUE)
empty_word <- readLines("character_string_for_empty_wells.txt")
threshold <- readLines("threshold_for_missingness_per_individual.txt")
threshold <- as.numeric(threshold)

data$content <- rep("DNA", nrow(data))
emptyLines <- grep(empty_word, data$INDV)
data[emptyLines, "content"] <- empty_word

#subset the data frame keeping only lines that content DNA (eclude empty wells)
data_noEmpty <- data[data$content == "DNA",]

#determine the outlier threshold
s <- summary(data_noEmpty$F_MISS)
q1 <- s[2]
q3 <- s[5]

IQR3 <- q3 + 3*(q3-q1)
IQR1.5 <- q3 + 1.5*(q3-q1)

outlier_threshold <- threshold

# tag the lines over the threshold
data$over_threshold <- rep("NO", nrow(data))
data[data$F_MISS >= outlier_threshold, "over_threshold"] <- "YES"
```

#### Distribution of missing genotypes for the identification of lines (taxa) to remove based on this criterion
###### Done at this filtering step:FastGBS_platypus_onlyPASS_noScaf_rmIndels_mac1_minMaxAll2_mac4_mm02.vcf
###### The filtering thresold of missing data was set manually to 0.8 (lines over 0.8 were removed)
<br/>
<br/>
```{r, echo = FALSE}
# plotting
ggplot(data = data, mapping = aes(y = F_MISS, x ="")) +
  geom_jitter(mapping = aes(colour = content, shape = over_threshold), height = 0) +
  geom_hline(aes(yintercept = c(outlier_threshold), linetype = paste("Manual threshold:", outlier_threshold)), colour = 'red') +
  geom_hline(aes(yintercept = IQR3, linetype = "Q3 + 3IQR"), colour = 'black') +
  geom_hline(aes(yintercept = IQR1.5, linetype = "Q3 + 1.5IQR"), colour = 'blue') +
  scale_linetype_manual(name = "", values = c(2,2,2), guide = guide_legend(override.aes = list(color = c("red", "blue", 'black')))) +
  ggtitle("Proportion of missing data per sample") +
  xlab("") + ylab("Proportion of missing data")
```
<br/>
<br/>

```{r, echo = FALSE}
to_remove <- data[data$over_threshold == "YES", "INDV"]
```
<br/>
<br/>

#### Printing the the barcodes/well/sampleName that were removed for being over the 0.8 missing data threshold
```{r, echo = FALSE}
print(to_remove, max.levels = 0)
```




```{r, echo = FALSE}
temp2 <- readLines("hetero_site_filt_input_name.txt")
data <- read.table(temp2, header = T)

matrix <- as.matrix(data[ ,3:ncol(data)])
matrix <- gsub("0/1", "1/0", matrix)
matrix <- matrix == "1/0"
matrix <- matrix*1
Proportion.Heterozygous <- rep(NA, nrow(matrix))

for (i in 1:nrow(matrix)) {
  Proportion.Heterozygous[i] <- sum(matrix[i,])/ncol(matrix)
}

data <- cbind(data, Proportion.Heterozygous)

#calculating the outliers based on the interquartile range
QT <- summary(data$Proportion.Heterozygous)[c(2, 5)]
QT <- as.vector(QT)
Q1 <- QT[1]
Q3 <- QT[2]
IQR3 <- Q3+(3*(Q3-Q1))
threshold <- Q3+(1.5*(Q3-Q1))

# filtering sites over Q3+1.5IQR
data_export <- data[data$Proportion.Heterozygous > threshold, c(1, 2)]
data <- data[, c("CHROM", "POS", "Proportion.Heterozygous")]

# producing the output file
write.table(data_export, file = "Heterozygosity_Outliers_sites_Q3plus1.5IQR.txt", col.names = F, row.names = F, quote = F)

# plotting
total_n_sites <- length(data$Proportion.Heterozygous)
over1.5IQR <- data$Proportion.Heterozygous >= threshold
o <- which(over1.5IQR)
n_sites_over_1.5IQR <- length(o)
over3IQR <- data$Proportion.Heterozygous >= IQR3
p <- which(over3IQR)
n_sites_over_3IQR <- length(p) 
n_section1 <- total_n_sites-n_sites_over_1.5IQR
perct_section1 <- (n_section1*100)/total_n_sites
p1 <- round(perct_section1, digits = 1)
n_section3 <- n_sites_over_3IQR
perct_section3 <- (n_section3*100)/total_n_sites
p3 <- round(perct_section3, digits = 1)
n_section2 <- total_n_sites-n_section1-n_section3
perct_section2 <- (n_section2*100)/total_n_sites
p2 <- round(perct_section2, digits = 1)
threshold <- as.numeric(format(round(threshold, 2), nsmall = 2))
IQR3 <- as.numeric(format(round(IQR3, 2), nsmall = 2))
nb_snps_tot <- length(data$Proportion.Heterozygous)
nb_snps_lost_1.5IQR <- length(which(data$Proportion.Heterozygous > threshold))
nb_snps_lost_3IQR <- length(which(data$Proportion.Heterozygous > IQR3))
nb_snps_lost_0.1 <- length(which(data$Proportion.Heterozygous > 0.1))

```
<br/>
<br/>
<br/>
<br/>

### Statistics made to filter out sites with high heterozygosity 
##### The threshold used for filtration is Q3+1.5IQR
##### The distribution was ploted at this filtering step:FastGBS_platypus_onlyPASS_noScaf_rmIndels_mac1_minMaxAll2_highNrm_mac4_mm20_imputed

<br/>
<br/>
```{r, echo = FALSE}
library(ggplot2)
ggplot(data = data, mapping = aes(y = Proportion.Heterozygous, x ="")) +
  geom_jitter(shape = 1, width = 1, color = "grey") +
  geom_violin(color = "black", alpha = 0.2) +
  geom_hline(yintercept = c(threshold, IQR3, 0.1), show.legend = TRUE, color = "black") +
  ggtitle("Proportion heterozygosity per SNP") +
  xlab("Sites (SNPs)") + ylab("Proportion heterozygosity") +
  geom_text(aes(0.1 , threshold, label = paste(threshold, "(Q3 + 1.5IQR)", ", filtering out", nb_snps_lost_1.5IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 0), color = "black") +
  geom_text(aes(0.1, IQR3, label = paste(IQR3,"(Q3 + 3IQR)", ", filtering out", nb_snps_lost_3IQR, "SNPs"), vjust = 0, hjust = -0.1, angle = 0), color = "black") +
  geom_text(aes(0.1, 0.1, label = paste("0.1", ", filtering out", nb_snps_lost_0.1, "SNPs"), vjust = 0, hjust = -0.3, angle = 0), color = "black") +
  geom_text(aes(0.1, 0.9, label = paste("Total nb of sites =", nb_snps_tot), vjust = 0, hjust = -2, angle = 0), color = "black")

```
<br/>
<br/>

###### sites over the threshold of Q3 + 1.5IQR were removed of the data set

###### The sites that were removed can be found in this output file: sites_Outliers_Heterozygosity_Q3plus1.5IQR.txt

<br/>
<br/>
<br/>
<br/>
<br/>
########################################################################################################################
# General statistics for the complete dataset after each filtration steps
<br/>

```{r, echo = FALSE}
overall_reports <- list.files(pattern = "Tassel_oall_report1")

# importing all Tassel overall reports to a unique data frame
overall_report_allvcfs <- NULL
for (i in overall_reports) {
  dat <- read.table(i, header=T, sep="\t", na.strings="", fill = F, stringsAsFactors = FALSE, dec = ",")
  name <- substr(i, 3, nchar(i)-24)
  colnames(dat)[2] <- name
  overall_report_allvcfs <- cbind(overall_report_allvcfs, dat[, name])
  colnames(overall_report_allvcfs)[ncol(overall_report_allvcfs)] <- name
  rownames(overall_report_allvcfs) <- dat$Stat.Type
}

# some data wrangling
overall_report_allvcfs <- as.data.frame(overall_report_allvcfs)
colnames(overall_report_allvcfs) <- gsub(" ", "_", colnames(overall_report_allvcfs))
rownames(overall_report_allvcfs) <- gsub(" ", "_", rownames(overall_report_allvcfs))
#colnames(overall_report_allvcfs)[3] <- "highNrm_taxa" 
overall_report_allvcfs <- t(overall_report_allvcfs)
d <- data.frame(rownames(overall_report_allvcfs))
overall_report_allvcfs <- cbind(overall_report_allvcfs, d)
overall_report_allvcfs <- overall_report_allvcfs[, c(ncol(overall_report_allvcfs), 1:(ncol(overall_report_allvcfs)-1))]
colnames(overall_report_allvcfs)[1] <- "filter"
filter_names_odered <- substr(overall_reports, 3, nchar(overall_reports)-24)
overall_report_allvcfs$filter <- factor(overall_report_allvcfs$filter, levels = filter_names_odered)
# dataframe type 2 for missing and non missing sites
df_number_not_missing <- overall_report_allvcfs[, c("filter", "Number_Not_Missing")]
df_number_not_missing$data_type <- rep("number_not_missing", nrow(df_number_not_missing))
colnames(df_number_not_missing)[2] <- "SNPs_number"
df_Number_Missing <- overall_report_allvcfs[, c("filter", "Number_Missing")]
df_Number_Missing$data_type <- rep("number_missing", nrow(df_Number_Missing))
colnames(df_Number_Missing)[2] <- "SNPs_number"
df_miss_not_miss <- rbind(df_Number_Missing, df_number_not_missing)
overall_report_allvcfs$Number_of_Sites <- as.numeric(overall_report_allvcfs$Number_of_Sites)
```

```{r, echo =FALSE}
knitr::kable(overall_report_allvcfs[, c(2, 3, 8, 15)], align = "cccc", caption = "")
```
<br/>
<br/>
<br/>
```{r, echo = FALSE}
# Number of SNPs at different filtering steps
ggplot(data = overall_report_allvcfs, mapping = aes(x = as.numeric(filter), y = Number_of_Sites)) +
  geom_point() +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  #scale_y_discrete(limits = rev) +
  labs(title = "Number of SNPs at different filtering steps", x = "Filters")
  
# plots for missing and non missing genotypes
#ggplot(data = df_miss_not_miss, mapping = aes(x = as.numeric(filter), y = SNPs_number, color = data_type)) +
 #geom_point() +
  #scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  #labs(title = "Stats on genotypes after different filtering steps", x = "Filters")
```
<br/>
<br/>
<br/>
```{r, echo =FALSE}
# plot for proportion missing data
ggplot(data = overall_report_allvcfs, mapping = aes(x = as.numeric(filter), y = Proportion_Missing)) +
  geom_point() +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion of missing data after different filtering steps", x = "Filters")
```
<br/>
<br/>
<br/>
```{r, echo =FALSE}
# plot for proportion heterozygous genotypes
ggplot(data = overall_report_allvcfs, mapping = aes(x = as.numeric(filter), y = Proportion_Heterozygous)) +
  geom_point() +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion of heterozygousity after different filtering steps", x = "Filters")
```

```{r, echo = FALSE}
############ importing the taxa summary by Tassel

taxa_reports <- list.files(pattern = "Tassel_taxa_report")

# importing all Tassel taxa reports to a unique data frame
taxa_report_allvcfs <- NULL
for (i in taxa_reports) {
  dat <- read.table(i, header=T, sep="\t", na.strings="", fill = F, stringsAsFactors = FALSE, dec = ",")
  name <- substr(i, 3, nchar(i)-24)
  dat$filters <- rep(name, nrow(dat))
  taxa_report_allvcfs <- rbind(taxa_report_allvcfs, dat)
}

# some data wrangling
taxa_report_allvcfs$filters <- as.factor(taxa_report_allvcfs$filters)
filter_names_odered <- substr(taxa_reports, 3, nchar(taxa_reports)-24)
taxa_report_allvcfs$filters <- factor(taxa_report_allvcfs$filters, levels = filter_names_odered)
taxa_report_allvcfs$Taxa <- gsub("^*", "s", taxa_report_allvcfs$Taxa)
```

<br/>
<br/>
<br/>
<br/>
```{r, echo= FALSE}
ggplot(data = taxa_report_allvcfs, mapping =aes(x = as.numeric(filters), y = as.numeric(Proportion.Missing), group = as.numeric(filters))) +
  geom_boxplot(outlier.shape = 3) + 
  geom_jitter(width = 0.1, alpha = 0) +
  geom_text(aes(label = Taxa), position = position_jitter(width=0.4,height=0.2), size = 1.5) +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion missing data per taxa after different filtering steps", x = "Filters") +
  scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
```
<br/>
#### Proportion missing values on the Y axis are expanded to help visualisation, real values can be found in this file: tassel_taxa_report_after_each_filtration_steps.txt
<br/>
<br/>
<br/>
<br/>
```{r, echo = FALSE}
ggplot(data = taxa_report_allvcfs, mapping =aes(x = as.numeric(filters), y = as.numeric(Proportion.Heterozygous), group = as.numeric(filters))) +
  geom_boxplot(outlier.shape = 3) + 
  geom_jitter(width = 0.1, alpha = 0) +
  geom_text(aes(label = Taxa), position = position_jitter(width=0.4,height=0.2), size = 1.5) +
  scale_x_continuous(breaks = 1:length(filter_names_odered), labels = filter_names_odered) +
  labs(title = "Proportion heterozygousity data per taxa after different filtering steps", x = "Filters") +
  scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
```
<br/>
#### Proportion heterozygousity values on the Y axis are expanded to help visualisation, real values can be found in this file: tassel_taxa_report_after_each_filtration_steps.txt
<br/>
<br/>
```{r, echo = FALSE}
dath <- taxa_report_allvcfs[as.character(taxa_report_allvcfs$filters) == "impute",]

ggplot(data = dath, mapping = aes(x = as.numeric(dath$Proportion.Heterozygous))) +
  geom_histogram(bins = 30) +
  labs(title = "Proportion heterozygousity per taxa in final imputed VCF file", x = "Proportion heterozygosity", y = "Count")
```
<br/>
<br/>
###Proportion heterozygousity statistics per taxa in final imputed VCF file
```{r, echo = FALSE}
print(summary(as.numeric(dath$Proportion.Heterozygous)))
```
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
#### The statistics and data points labels correspondance can be found in the following output file:tassel_taxa_report_after_each_filtration_steps.txt

' > Final_report_all_stats_graphs.Rmd

# knitting the rmarkdown script to output the html graphs report
R -e "rmarkdown::render('Final_report_all_stats_graphs.Rmd')"


# removing temporary files and scripts
#rm Final_report_all_stats_graphs.Rmd
rm Rplots.pdf
#rm Rscript_tassel_overall_and_taxa_reports_for_all_filtration_steps.R
#rm Rscript_imiss_stats_to_remove_lines_over_missing_data_threshold.R
#rm temp_input_name_for_R.txt
#rm hetero_site_filt_input_name.txt
#rm *report.txt
#rm *report1.txt
#rm *report2.txt
