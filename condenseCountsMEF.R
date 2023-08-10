setwd("~/Documents/GradSchool/WekLab/MikeData/tRNA-Seq")
library(readr)
library(dplyr)

#####create empty data frame to hold counts
#retreive tRNA names from mm10-tRNAsPheCCA.fa file
my_genes_file <- read_lines("mm10-tRNAsCCA.fa")
#take only files containing reference gene names
my_genes_pre <- my_genes_file[seq(1, 452, 2)]
#declare empty vector to hold values
my_gene_names <- c()
#iterate over 
for (x in 1:227) {
     split_temp <- strsplit(as.character(my_genes_pre[x]), split = "-", fixed = TRUE)
     name_temp <- paste(split_temp[[1]][2], split_temp[[1]][3], sep = "-")
     my_gene_names <- c(my_gene_names, name_temp)
}
#Remove duplicate entries
unique_gene_names <- unique(my_gene_names)
#create empty list of data frames to hold counts
my_counts <- data.frame(gene_names=unique_gene_names, total=rep(0,52), charged=rep(0,52), uncharged=rep(0,52), chargedPercent=rep(0,52))
rownames(my_counts) <- unique_gene_names
my_countFiles <- c("cont1.counts", "cont2.counts", "cont3.counts", "HF1.counts", "HF2.counts", "HF3.counts")
isodecoder_counts <- list()
for (k in 1:6){
     isodecoder_counts[[my_countFiles[k]]] <- my_counts
}

####read counts table produced by parseSamCharge.py into table. Must set this path to the directory containing the counts files, here the director is "FullCounts"
setwd("~/Documents/GradSchool/WekLab/MikeData/tRNA-Seq/2021Apr28JagannathData/FullCounts")
count_tables <- list()
for (filename in my_countFiles){
     my_data<- read.table(filename, header = FALSE, sep = ":", strip.white = TRUE)
     my_data2 <- cbind(my_data, rep(0,227))
     my_tbl <- as_tibble(my_data2)
     colnames(my_tbl) <- c("gene_names", "total", "charged", "uncharged", "chargedPercent")
     my_tbl <- my_tbl %>% mutate(chargedPercent = charged/total)
     my_tbl_sorted <- arrange(my_tbl, gene_names)
     count_tables[[as.character(filename)]] <- my_tbl_sorted
}
     



##combine counts of different genes of same isodecoder
for (j in seq(1,6)){
     for (i in 1:length(count_tables[[j]]$gene_names)){
          gene_name_tmp <- as.character(count_tables[[j]]$gene_names[i])
          name_split_tmp <- strsplit(gene_name_tmp, split = "-", fixed = TRUE)
          new_name_tmp <- paste(name_split_tmp[[1]][2], name_split_tmp[[1]][3], sep = "-")
          if (new_name_tmp %in% isodecoder_counts[[j]][, 1]) {
               isodecoder_counts[[j]][new_name_tmp, 2] = (isodecoder_counts[[j]][new_name_tmp, 2] + count_tables[[j]][i, 2])
               isodecoder_counts[[j]][new_name_tmp, 3] = (isodecoder_counts[[j]][new_name_tmp, 3] + count_tables[[j]][i, 3])
               isodecoder_counts[[j]][new_name_tmp, 4] = (isodecoder_counts[[j]][new_name_tmp, 4] + count_tables[[j]][i, 4])
          } 
          isodecoder_counts[[j]][,5] <- isodecoder_counts[[j]][,3]/(isodecoder_counts[[j]][,3]+isodecoder_counts[[j]][,4])
     }
}
#normalize data
total_counts <- c()
for (k in seq(1,6)){
     chargeAndUncharge <- sum(count_tables[[k]]$charged[2:227])+sum(count_tables[[k]]$uncharged[2:227])
     total_counts <- c(total_counts, chargeAndUncharge)
}
average_count <- sum(total_counts)/6
norm_const <- total_counts/average_count
#create list of dataframes with normalized data
normalized_isodecoder_counts <- list()
for (m in seq(1,6)){
     names <- unique_gene_names
     cc <- isodecoder_counts[[m]]$charged/norm_const[m]
     uc <- isodecoder_counts[[m]]$uncharged/norm_const[m]
     pc <- cc/(cc+uc)
     normalized_table <- data.frame(gene_names=names,Charged_Count=cc,Uncharged_Count=uc,Percent_Charged=pc)
     colnames(normalized_table) <- c("gene_names", "Charged_Count", "Uncharged_Count", "Percent_Charged")
     normalized_isodecoder_counts[[my_countFiles[m]]] <- normalized_table
}

#create list to hold averages
norm_avg <- list()
norm_cont <- data.frame(gene_names=unique_gene_names, cont1_charged=normalized_isodecoder_counts[[1]]$Charged_Count, cont2_charged=normalized_isodecoder_counts[[2]]$Charged_Count, cont3_charged=normalized_isodecoder_counts[[3]]$Charged_Count,
                        cont1_uncharged=normalized_isodecoder_counts[[1]]$Uncharged_Count, cont2_uncharged=normalized_isodecoder_counts[[2]]$Uncharged_Count, cont3_uncharged=normalized_isodecoder_counts[[3]]$Uncharged_Count,
                        cont_charged_avg=rep(0,52), cont_charged_sd=rep(0,52), cont_uncharged_avg=rep(0,52), cont_uncharged_sd=rep(0,52), cont_chargePercent_avg=rep(0,52), cont_chargePercent_sd=rep(0,52), cont_totalSD=rep(0,52))
norm_avg[["control"]] <- norm_cont
norm_HF <- data.frame(gene_names=unique_gene_names, HF1_charged=normalized_isodecoder_counts[[4]]$Charged_Count, HF2_charged=normalized_isodecoder_counts[[5]]$Charged_Count, HF3_charged=normalized_isodecoder_counts[[6]]$Charged_Count,
                      HF1_uncharged=normalized_isodecoder_counts[[4]]$Uncharged_Count, HF2_uncharged=normalized_isodecoder_counts[[5]]$Uncharged_Count, HF3_uncharged=normalized_isodecoder_counts[[6]]$Uncharged_Count,
                      HF_charged_avg=rep(0,52), HF_charged_SD=rep(0,52), HF_uncharged_avg=rep(0,52), HF_uncharged_SD=rep(0,52), HF_chargePercent_avg=rep(0,52), HF_chargePercent_sd=rep(0,52), HF_totalSD=rep(0,52))
norm_avg[["HF"]] <- norm_HF

#calculate average and sd on normalized data
library(genefilter)
sample1_charge <- c()
sample2_charge <- c()
sample3_charge <- c()
total1 <- c()
total2 <- c()
total3 <- c()
chargeSD <- function(dex){
     return(sd(c(sample1_charge[dex], sample2_charge[dex], sample3_charge[dex])))    
}
chargeAvg <- function(dex){
     return((sample1_charge[dex]+sample2_charge[dex]+sample3_charge[dex])/3)
}
totalSD <- function(dex){
     return(sd(c(total1[dex], total2[dex], total3[dex])))
}
for (n in seq(1,2)){
     norm_avg[[n]][,8] <- (norm_avg[[n]][,2] + norm_avg[[n]][,3] + norm_avg[[n]][,4])/3
     norm_avg[[n]][,9] <- rowSds(as.matrix(norm_avg[[n]][,c(2,3,4)]))
     norm_avg[[n]][,10] <- (norm_avg[[n]][,5] + norm_avg[[n]][,6] + norm_avg[[n]][,7])/3
     norm_avg[[n]][,11] <- rowSds(as.matrix(norm_avg[[n]][,c(5,6,7)]))
     norm_avg[[n]][,12] <- norm_avg[[n]][,8]/(norm_avg[[n]][,8]+norm_avg[[n]][,10])
     
     sample1_charge <- norm_avg[[n]][,2]/(norm_avg[[n]][,2]+norm_avg[[n]][,5])
     sample2_charge <- norm_avg[[n]][,3]/(norm_avg[[n]][,3]+norm_avg[[n]][,6])
     sample3_charge <- norm_avg[[n]][,4]/(norm_avg[[n]][,4]+norm_avg[[n]][,7])
     
     total1 <- norm_avg[[n]][,2]+norm_avg[[n]][,5]
     total2 <- norm_avg[[n]][,3]+norm_avg[[n]][,6]
     total3 <- norm_avg[[n]][,4]+norm_avg[[n]][,7]
     
     for (j in 1:length(sample1_charge)){
          norm_avg[[n]][j, 12] <- chargeAvg(j)
          norm_avg[[n]][j, 13] <- chargeSD(j)
          norm_avg[[n]][j, 14] <- totalSD(j)
     }
}

#remove NaN values from norm_avg tables
norm_avg2 <- list()
norm_avg2[["control"]] <- norm_avg[["control"]][-c(42,52),]
norm_avg2[["HF"]] <- norm_avg[["HF"]][-c(42,52),]


#remove rows with NaN values for charge percent SD also change labels of control and starvation samples
norm_avg3 <- list()
norm_avg3[["control"]] <- norm_avg2[["control"]][-c(39,41), ]
norm_avg3[["HF"]] <- norm_avg2[["HF"]][-c(39,41), ]

#prepare data for heatmap
library(ggplot2)

#plot control and HF only
y <- as.character(norm_avg3[[1]][,1]) #this is a character vector of the gene names
x <- c("Control", "HF")
data = expand.grid(X=x,Y=y)
data$charge_percent <- as.vector(rbind(norm_avg3[[1]][,12], norm_avg3[[2]][,12]))
data$charge_percentSD <- as.vector(rbind(norm_avg3[[1]][,13], norm_avg3[[2]][,13]))
colnames(data) <- c("treatment", "gene_name", "charge_percent", "charge_percentSD")
heatmap <- ggplot(data, aes(x=treatment,y=factor(gene_name, levels = rev(levels(factor(gene_name)))),fill = charge_percent))+geom_tile() + scale_fill_gradient(low = "red", high = "blue") +
     labs(title = "tRNA Charging", x="Treatment", y="Isodecoder", fill = "Fraction\nCharged") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0)) + coord_equal() + scale_x_discrete(labels = c("Control", "HF"))
#ggsave("HeatmapIsodecoderHF.png", width = 1200, height = 2400, units = "px", dpi = 300)

my_charge_barplot <- ggplot(data, aes(fill=treatment, y=charge_percent, x=gene_name, group = treatment)) + geom_bar(position = "dodge", stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     geom_errorbar(aes(ymin=charge_percent-charge_percentSD, ymax=charge_percent+charge_percentSD), width = 0.2, position = position_dodge(0.9)) +
     labs(title = "tRNA Charging", y="Fraction Charged", x="Isodecoder")
#ggsave("IsodecoderChargingHF.png", plot = my_charge_barplot, width = 2400, height = 1200, units = "px", dpi = 300)

# calculate statistics
source('~/Documents/GradSchool/WekLab/t.test.Welch.R')
HF_pvalues <- c()

for (j in 1:length(norm_avg3[[1]]$gene_names)){
     HF_pvalues <- c(HF_pvalues, t.test.Welch(norm_avg3[[1]]$cont_chargePercent_avg[j], norm_avg3[[2]]$HF_chargePercent_avg[j], norm_avg3[[1]]$cont_chargePercent_sd[j], norm_avg3[[2]]$HF_chargePercent_sd[j], 3, 3)[5])
}
cont_pvalues <- rep(1,48)
data$pvalues <- as.vector(rbind(cont_pvalues, HF_pvalues))
adj_pvalues <- c()
for (i in 1:length(data$pvalues)){
     adj_tmp <- data$pvalues[i]*48
     if (adj_tmp < 1){
          adj_pvalues <- c(adj_pvalues, adj_tmp)
     } else {
          adj_pvalues <- c(adj_pvalues, 1.000000000)
     }
}
data$padj <- adj_pvalues

#calculate adjusted P values
#HF 
HF_stats <- data.frame(gene_names=norm_avg3[["HF"]]$gene_names, chargePercent=norm_avg3[["HF"]]$HF_chargePercent_avg, p_value=HF_pvalues)
HF_stats_tbl <- as_tibble(HF_stats)
HF_stats_tbl_sorted <- arrange(HF_stats_tbl, p_value)
HF_stats_tbl_sorted$rank <- seq(1:48)
HF_stats_tbl_sorted$BHFDR <- HF_stats_tbl_sorted$p_value*48/HF_stats_tbl_sorted$rank



#######################################
#combined amino acid frequencies
#######################################
#load user defined t.test.Welch function
source('~/Documents/GradSchool/WekLab/t.test.Welch.R')
codon_table <- list()
codon_table[["ile"]] <- c("ATT", "ATC", "ATA")
codon_table[["leu"]] <- c("CTT","CTG", "CTC", "CTA", "TTA", "TTG")
codon_table[["val"]] <- c("GTT", "GTC", "GTA", "GTG")
codon_table[["phe"]] <- c("TTT", "TTC")
codon_table[["met"]] <- c("ATG")
codon_table[["cys"]] <- c("TGT", "TGC")
codon_table[["ala"]] <- c("GCA", "GCC", "GCG", "GCT")
codon_table[["gly"]] <- c("GGA", "GGC", "GGG", "GGT")
codon_table[["pro"]] <- c("CCA", "CCC", "CCG", "CCT")
codon_table[["thr"]] <- c("ACA", "ACC", "ACG", "ACT")
codon_table[["ser"]] <- c("TCA", "TCC", "TCG", "TCT", "AGT", "AGC")
codon_table[["tyr"]] <- c("TAT", "TAC")
codon_table[["trp"]] <- c("TGG")
codon_table[["gln"]] <- c("CAA", "CAG")
codon_table[["asn"]] <- c("AAT", "AAC")
codon_table[["his"]] <- c("CAC", "CAT")
codon_table[["glu"]] <- c("GAA", "GAG")
codon_table[["asp"]] <- c("GAC", "GAT")
codon_table[["lys"]] <- c("AAA", "AAG")
codon_table[["arg"]] <- c("CGA", "CGC", "CGG", "CGT", "AGA", "AGG")
codon_table[["stop"]] <- c("TAA", "TAG", "TGA")

#make data frame to hold values for p-site amino acids
aa_stats <- as.data.frame(matrix(NA, ncol = 6, nrow = 21))
colnames(aa_stats) <- c("amino_acid", "control_avg", "control_sd", "hf_avg", "hf_sd","hf_p.value")
rownames(aa_stats) <- names(codon_table)
#loop over table to calculate Average for each codon then sum them and calculate the stdev
for (i in 1:length(codon_table)){
     #calculate averages
     control_vect <- norm_avg3[[1]][grep(pattern = names(codon_table)[i], x = norm_avg3[[1]]$gene_names, ignore.case = TRUE), 12]
     hf_vect <- norm_avg3[[2]][grep(pattern = names(codon_table)[i], x = norm_avg3[[1]]$gene_names, ignore.case = TRUE), 12]
     

     control_mean <- sum(control_vect)/length(control_vect)
     hf_mean <- sum(hf_vect)/length(hf_vect)
     
     
     #calculate SD
     control_vectsd <- norm_avg3[[1]][grep(pattern = names(codon_table)[i], x = norm_avg3[[1]]$gene_names, ignore.case = TRUE), 13]
     hf_vectsd <- norm_avg3[[2]][grep(pattern = names(codon_table)[i], x = norm_avg3[[1]]$gene_names, ignore.case = TRUE), 13]
    
     control_sd <- sqrt(sum(control_vectsd^2))/length(control_vectsd)
     hf_sd <- sqrt(sum(hf_vectsd^2))/length(hf_vectsd)
  
     
     #add values to data frame
     aa_stats[i, 2] <- control_mean
     aa_stats[i, 3] <- control_sd
     aa_stats[i, 4] <- hf_mean
     aa_stats[i, 5] <- hf_sd
  
     pval_temp <- t.test.Welch(m1=control_mean, m2=hf_mean, s1=control_sd, s2=hf_sd, 3,3)
     aa_stats[i, 6] <- pval_temp["p.value"]

}
aa_stats$amino_acid <- rownames(aa_stats)

#calculate adjusted P values
#HF treatment
aa_stats_tbl <- as_tibble(aa_stats[1:20, ])
HF_aa_stats_tbl_sorted <- arrange(aa_stats_tbl, hf_p.value)
HF_aa_stats_tbl_sorted$HF_rank <- seq(1:20)
HF_aa_stats_tbl_sorted$HF_BHFDR <- HF_aa_stats_tbl_sorted$hf_p.value*20/HF_aa_stats_tbl_sorted$HF_rank


#Arrange results table
aa_stats_tbl <- arrange(aa_stats_tbl, amino_acid)


#Plot control and HF
y = aa_stats_tbl$amino_acid
x = c("Control", "HF")
data = expand.grid(X=x,Y=y)
data$charge_percent <- as.vector(rbind(aa_stats_tbl$control_avg, aa_stats_tbl$hf_avg))
data$charge_percentSD <- as.vector(rbind(aa_stats_tbl$control_sd, aa_stats_tbl$hf_sd))
data$sig <- c(rep("",29), "*",rep("",10))
colnames(data) <- c("treatment", "amino_acid", "charge_percent", "charge_percentSD", "sig")

heatmap <- ggplot(data, aes(x=treatment,y=factor(amino_acid, levels = rev(levels(factor(amino_acid)))),fill = charge_percent))+geom_tile() + scale_fill_gradient(low = "red", high = "blue") +
     labs(title = "tRNA Charging", x="Treatment", y="Amino Acid", fill = "Fraction\nCharged") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0)) + coord_equal() + scale_x_discrete(labels = c("Control", "HF"))
#ggsave("HeatmapAminoAcidHF.png", width = 1200, height = 2400, units = "px", dpi = 300)

my_charge_barplot <- ggplot(data, aes(fill=treatment, y=charge_percent, x=amino_acid, group = treatment)) + geom_bar(position = "dodge", stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     geom_errorbar(aes(ymin=charge_percent-charge_percentSD, ymax=charge_percent+charge_percentSD), width = 0.2, position = position_dodge(0.9)) +
     labs(title = "tRNA Charging", y="Fraction Charged", x="Amino Acid")

#To save delete comment marker below:
#ggsave("AminoAcidChargingHF.png", plot = my_charge_barplot, width = 2400, height = 1200, units = "px", dpi = 300)


