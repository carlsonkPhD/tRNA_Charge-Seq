setwd("~/Documents/GradSchool/WekLab/MikeData/tRNA-Seq")
library(readr)
library(dplyr)

#####create empty data frame to hold counts
#retreive tRNA names from hg38tRNAsCCA.fa file
my_genes_file <- read_lines("hg38tRNAsCCA.fa")
#take only files containing reference gene names
my_genes_pre <- my_genes_file[seq(1, 520, 2)]
#declare empty vector to hold values
my_gene_names <- c()
#iterate over 
for (x in 1:260) {
     split_temp <- strsplit(as.character(my_genes_pre[x]), split = "-", fixed = TRUE)
     name_temp <- paste(split_temp[[1]][2], split_temp[[1]][3], sep = "-")
     my_gene_names <- c(my_gene_names, name_temp)
}
#Remove duplicate entries
unique_gene_names <- unique(my_gene_names)
#create empty list of data frames to hold counts
my_counts <- data.frame(gene_names=unique_gene_names, total=rep(0,49), charged=rep(0,49), uncharged=rep(0,49), chargedPercent=rep(0,49))
rownames(my_counts) <- unique_gene_names


####read counts table produced by parseSamCharge.py into table
setwd("~/Documents/GradSchool/WekLab/MikeData/tRNA-Seq/2021Aug4RicardoData/FullCountsHuman")
my_countFiles <- list.files()
isodecoder_counts <- list()
for (k in 1:length(my_countFiles)){
     isodecoder_counts[[my_countFiles[k]]] <- my_counts
}

count_tables <- list()
for (filename in my_countFiles){
     my_data<- read.table(filename, header = FALSE, sep = ":", strip.white = TRUE)
     my_data2 <- cbind(my_data, rep(0,261))
     my_tbl <- as_tibble(my_data2)
     colnames(my_tbl) <- c("gene_names", "total", "charged", "uncharged", "chargedPercent")
     my_tbl <- my_tbl %>% mutate(chargedPercent = charged/total)
     my_tbl_sorted <- arrange(my_tbl, gene_names)
     count_tables[[as.character(filename)]] <- my_tbl_sorted
}
     



##combine counts of different genes of same isodecoder
for (j in 1:length(my_countFiles)){
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
for (k in 1:length(my_countFiles)){
     chargeAndUncharge <- sum(count_tables[[k]]$charged[2:227])+sum(count_tables[[k]]$uncharged[2:227])
     total_counts <- c(total_counts, chargeAndUncharge)
}
average_count <- sum(total_counts)/length(my_countFiles)
norm_const <- total_counts/average_count
#create list of dataframes with normalized data
normalized_isodecoder_counts <- list()
for (m in 1:length(my_countFiles)){
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
norm_cont <- data.frame(gene_names=unique_gene_names, cont1_charged=normalized_isodecoder_counts[[1]]$Charged_Count, cont2_charged=normalized_isodecoder_counts[[2]]$Charged_Count, cont3_charged=normalized_isodecoder_counts[[3]]$Charged_Count, cont4_charged=normalized_isodecoder_counts[[4]]$Charged_Count,
                        cont1_uncharged=normalized_isodecoder_counts[[1]]$Uncharged_Count, cont2_uncharged=normalized_isodecoder_counts[[2]]$Uncharged_Count, cont3_uncharged=normalized_isodecoder_counts[[3]]$Uncharged_Count, cont4_uncharged=normalized_isodecoder_counts[[4]]$Uncharged_Count,
                        cont_charged_avg=rep(0,49), cont_charged_sd=rep(0,49), cont_uncharged_avg=rep(0,49), cont_uncharged_sd=rep(0,49), cont_chargePercent_avg=rep(0,49), cont_chargePercent_sd=rep(0,49), cont_totalSD=rep(0,49))
norm_avg[["control"]] <- norm_cont
norm_IB <- data.frame(gene_names=unique_gene_names, IB1_charged=normalized_isodecoder_counts[[5]]$Charged_Count, IB2_charged=normalized_isodecoder_counts[[6]]$Charged_Count, IB3_charged=normalized_isodecoder_counts[[7]]$Charged_Count, IB4_charged=normalized_isodecoder_counts[[8]]$Charged_Count,
                      IB1_uncharged=normalized_isodecoder_counts[[5]]$Uncharged_Count, IB2_uncharged=normalized_isodecoder_counts[[6]]$Uncharged_Count, IB3_uncharged=normalized_isodecoder_counts[[7]]$Uncharged_Count, IB4_uncharged=normalized_isodecoder_counts[[8]]$Uncharged_Count,
                      IB_charged_avg=rep(0,49), IB_charged_SD=rep(0,49), IB_uncharged_avg=rep(0,49), IB_uncharged_SD=rep(0,49), IB_chargePercent_avg=rep(0,49), IB_chargePercent_sd=rep(0,49), IB_totalSD=rep(0,49))
norm_avg[["IB"]] <- norm_IB
norm_IBAA <- data.frame(gene_names=unique_gene_names, IBAA_charged=normalized_isodecoder_counts[[9]]$Charged_Count, IBAA2_charged=normalized_isodecoder_counts[[10]]$Charged_Count, IBAA3_charged=normalized_isodecoder_counts[[11]]$Charged_Count, IBAA4_charged=normalized_isodecoder_counts[[12]]$Charged_Count,
                        IBAA1_uncharged=normalized_isodecoder_counts[[9]]$Uncharged_Count, IBAA2_uncharged=normalized_isodecoder_counts[[10]]$Uncharged_Count, IBAA3_uncharged=normalized_isodecoder_counts[[11]]$Uncharged_Count, IBAA4_uncharged=normalized_isodecoder_counts[[12]]$Uncharged_Count,
                        IBAA_charged_avg=rep(0,49), IBAA_charged_sd=rep(0,49), IBAA_uncharged_avg=rep(0,49), IBAA_uncharged_sd=rep(0,49), IBAA_chargePercent_avg=rep(0,49), IBAA_chargePercent_sd=rep(0,49), IBAA_totalSD=rep(0,49))
norm_avg[["IBAA"]] <- norm_IBAA
#calculate average and sd on normalized data
library(genefilter)
sample1_charge <- c()
sample2_charge <- c()
sample3_charge <- c()
sample4_charge <- c()
total1 <- c()
total2 <- c()
total3 <- c()
total4 <- c()
chargeSD <- function(dex){
     return(sd(c(sample1_charge[dex], sample2_charge[dex], sample3_charge[dex], sample4_charge[dex])))    
}
chargeAvg <- function(dex){
     return((sample1_charge[dex]+sample2_charge[dex]+sample3_charge[dex]+sample4_charge[dex])/4)
}
totalSD <- function(dex){
     return(sd(c(total1[dex], total2[dex], total3[dex], total4[dex])))
}
for (n in seq(1,3)){
     norm_avg[[n]][,10] <- (norm_avg[[n]][,2] + norm_avg[[n]][,3] + norm_avg[[n]][,4] + norm_avg[[n]][,5])/4 #charged_avg
     norm_avg[[n]][,11] <- rowSds(as.matrix(norm_avg[[n]][,c(2,3,4,5)])) #charged_SD
     norm_avg[[n]][,12] <- (norm_avg[[n]][,6] + norm_avg[[n]][,7] + norm_avg[[n]][,8] + norm_avg[[n]][,9])/4 #uncharged_avg
     norm_avg[[n]][,13] <- rowSds(as.matrix(norm_avg[[n]][,c(6,7,8,9)])) #uncharged_SD
     
     sample1_charge <- norm_avg[[n]][,2]/(norm_avg[[n]][,2]+norm_avg[[n]][,6])
     sample2_charge <- norm_avg[[n]][,3]/(norm_avg[[n]][,3]+norm_avg[[n]][,7])
     sample3_charge <- norm_avg[[n]][,4]/(norm_avg[[n]][,4]+norm_avg[[n]][,8])
     sample4_charge <- norm_avg[[n]][,5]/(norm_avg[[n]][,5]+norm_avg[[n]][,9])
     
     total1 <- norm_avg[[n]][,2]+norm_avg[[n]][,6]
     total2 <- norm_avg[[n]][,3]+norm_avg[[n]][,7]
     total3 <- norm_avg[[n]][,4]+norm_avg[[n]][,8]
     total4 <- norm_avg[[n]][,5]+norm_avg[[n]][,9] 
     
     for (j in 1:length(sample1_charge)){
          norm_avg[[n]][j, 14] <- chargeAvg(j)
          norm_avg[[n]][j, 15] <- chargeSD(j)
          norm_avg[[n]][j, 16] <- totalSD(j)
     }
}

#remove NaN values from norm_avg tables (make union of NaN values from all tables then remove those rows from all)
norm_avg2 <- list()
na1 <- which(is.na(norm_avg[[1]]$cont_chargePercent_avg))
na2 <- which(is.na(norm_avg[[2]]$IB_chargePercent_avg))
na3 <- which(is.na(norm_avg[[3]]$IBAA_chargePercent_avg))
naUnion <- unique(c(na1,na2,na3))
norm_avg2[["control"]] <- norm_avg[["control"]][-naUnion,]
norm_avg2[["IB"]] <- norm_avg[["IB"]][-naUnion,]
norm_avg2[["IBAA"]] <- norm_avg[["IBAA"]][-naUnion,]


#prepare data for heatmap
library(ggplot2)
y <- as.character(norm_avg2[[1]][,1])
x <- c("cont.", "IB", "IB+AA")
data <- expand.grid(X=x, Y=y)
data$charge_percent <- as.vector(rbind(norm_avg2[[1]][,14], norm_avg2[[2]][,14], norm_avg2[[3]][,14]))
data$charge_percentSD <- as.vector(rbind(norm_avg2[[1]][,15], norm_avg2[[2]][,15], norm_avg2[[3]][,15]))
data$total <- as.vector(rbind(norm_avg2[[1]][,10]+norm_avg2[[1]][,12], norm_avg2[[2]][,10]+norm_avg2[[2]][,12], norm_avg2[[3]][,10]+norm_avg2[[3]][,12]))
data$totalSD <- as.vector(rbind(norm_avg2[[1]][,16], norm_avg2[[2]][,16], norm_avg2[[3]][,16]))
colnames(data) <- c("treatment", "gene_name", "charge_percent", "charge_percentSD", "total", "totalSD")
heatmap <- ggplot(data, aes(treatment,gene_name, fill=charge_percent))+geom_tile() + scale_fill_gradient(low = "red", high = "blue") + labs(title = "Isodecoder Charge Percent", x="Percent Charged", y="Isodecoder")
my_charge_barplot <- ggplot(data, aes(fill=treatment, y=charge_percent, x=gene_name)) + geom_bar(position = "dodge", stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    geom_errorbar(aes(ymin=charge_percent-charge_percentSD, ymax=charge_percent+charge_percentSD), width = 0.2, position = position_dodge(0.9)) +
                    labs(title = "Isodecoder Charge Percent", y="Percent Charged", x="Isodecoder")
my_total_barplot <- ggplot(data, aes(fill=treatment, y=total, x=gene_name)) + geom_bar(position = "dodge", stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                    geom_errorbar(aes(ymin=total-totalSD, ymax=total+totalSD), width = 0.2, position = position_dodge(0.9)) +
                    labs(title = "Isodecoder Counts", y="Counts", x="Isodecoder")
                    
# calculate statistics
source('~/Documents/GradSchool/WekLab/t.test.Welch.R')
IB_pvalues <- c()
IBAA_pvalues <- c()
for (j in 1:length(norm_avg2[[1]]$gene_names)){
     IB_pvalues <- c(IB_pvalues, t.test.Welch(norm_avg2[[1]]$cont_chargePercent_avg[j], norm_avg2[[2]]$IB_chargePercent_avg[j], norm_avg2[[1]]$cont_chargePercent_sd[j], norm_avg2[[2]]$IB_chargePercent_sd[j], 4, 4)[5])
     IBAA_pvalues <- c(IBAA_pvalues, t.test.Welch(norm_avg2[[1]]$cont_chargePercent_avg[j], norm_avg2[[3]]$IBAA_chargePercent_avg[j], norm_avg2[[1]]$cont_chargePercent_sd[j], norm_avg2[[3]]$IBAA_chargePercent_sd[j], 4, 4)[5])
}
cont_pvalues <- rep(1,47)
data$pvalues <- as.vector(rbind(cont_pvalues, IB_pvalues, IBAA_pvalues))
adj_pvalues <- c()
for (i in 1:length(data$pvalues)){
     adj_tmp <- data$pvalues[i]*47
     if (adj_tmp < 1){
          adj_pvalues <- c(adj_pvalues, adj_tmp)
     } else {
          adj_pvalues <- c(adj_pvalues, 1.000000000)
     }
}
data$padj <- adj_pvalues

#calculate adjusted P values
#IB 
IB_stats <- data.frame(gene_names=norm_avg2[["IB"]]$gene_names, chargePercent=norm_avg2[["IB"]]$IB_chargePercent_avg, p_value=IB_pvalues)
IB_stats_tbl <- as_tibble(IB_stats)
IB_stats_tbl_sorted <- arrange(IB_stats_tbl, p_value)
IB_stats_tbl_sorted$rank <- seq(1:47)
IB_stats_tbl_sorted$BHFDR <- IB_stats_tbl_sorted$p_value*47/IB_stats_tbl_sorted$rank

IBAA_stats <- data.frame(gene_names=norm_avg2[["IBAA"]]$gene_names, chargePercent=norm_avg2[["IBAA"]]$IBAA_chargePercent_avg, p_value=IBAA_pvalues)
IBAA_stats_tbl <- as_tibble(IBAA_stats)
IBAA_stats_tbl_sorted <- arrange(IBAA_stats_tbl, p_value)
IBAA_stats_tbl_sorted$rank <- seq(1:47)
IBAA_stats_tbl_sorted$BHFDR <- IBAA_stats_tbl_sorted$p_value*47/IBAA_stats_tbl_sorted$rank

###correlation analysis
#construct data table where each column is a diff sample and each row is a tRNA gene
L <- list()
for (i in 1:length(my_countFiles)){
     L[[i]] <- isodecoder_counts[[i]][,5]
}
m <- do.call(cbind, L)
colnames(m) <- my_countFiles
rownames(m) <- isodecoder_counts[[1]][,1]
new_m <- m[-naUnion, ]
res <- cor(new_m)
round_res <- round(res, 2)
