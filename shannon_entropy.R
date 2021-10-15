#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args[1] = input_path # input/path/to/MSA_file/of/an/HOG
# args[2]= out_path # output/path
library(msa)
library(entropy)
library(tidyverse)
library(ggplot2)
library(cowplot)

hvSiteEntCutoff <- 1.5
MinGapFraction <- 1
MinGapBlockWidth <- 1
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

file_name=args[1]
total_msa <- readAAMultipleAlignment(file_name)

HOG <- str_split(file_name, "\\.")[[1]][7]
print(paste("Looking at HOG",HOG,"..."))

gene_list <- rownames(total_msa)

# seperate sequences as status
cult_list <- c("GLA4","GP3","GP22","GP51","GP72",
               "GP772-1","GP295-1","GP39","GP77","GP536",
               "GP640","GP761-1","Nipponbare","GP551","GP567",
               "GP669","GP677","WYG7","KY131","DHX2",
               "IL9","Koshihikari","LG31")

land_list <- c("HP119","HP263","HP274","HP327","HP362-2",
               "HP383","HP396","HP407","HP486","HP492",
               "HP517-1","HP577","GP104","GP124","GP540",
               "GP62","Kasalath","HP13-2","HP14","HP38",
               "HP44","HP45","HP48","HP91-2","HP98",
               "HP103","HP314","HP390","UR28")

wild_list <- c("W0123-1","W0141","W0170","W1687","W1698",
               "W1739","W1754","W1777","W1943","W1979",
               "W2012","W3078-2","W3095-2")
feral_list <- c("W0128","W3105-1")

cult_sub <- c()
land_sub <- c()
wild_sub <- c()
feral_sub <- c()

for (i in gene_list) {
    out <- gsub("_.*", "", i)
    index <- which(gene_list == i)
    if (out %in% cult_list) {
        cult_sub <- append(cult_sub, index)
        }
    else if (out %in% land_list) {
        land_sub <- append(land_sub, index)
        }
    else if (out %in% wild_list) {
        wild_sub <- append(wild_sub, index)
        }
    else if (out %in% feral_list) {
        feral_sub <- append(feral_sub, index)
    }
    else {NULL}
}

# subset MSA according to status
cult_msa <- total_msa
rowmask(cult_msa) <- IRanges(c(land_sub, wild_sub, feral_sub))
land_msa <- total_msa
rowmask(land_msa) <- IRanges(c(cult_sub, wild_sub, feral_sub))
wild_msa <- total_msa
rowmask(wild_msa) <- IRanges(c(cult_sub, land_sub, feral_sub))

cal_entropy <- function(MSA, MGF, MGBW, ALPHA_21, ALPHA_20, HOG_ID, HSEC){
#     if ("-" %in% rownames(consensusMatrix(MSA))){
#     autoMasked <- maskGaps(MSA)
# #                            min.fraction = MGF,
# #                            min.block.width = MGBW) ##KEY FILTERING PARAMETERS
#     MinAli <- as(autoMasked, "AAStringSet")
#     }else{
    MinAli<-as(MSA, "AAStringSet")
    
     ## Calculating Consensus Matrix
    Tidy_CM<-as_tibble(t(consensusMatrix(MinAli, baseOnly = T)))
    
    ## Compensating for consensus matrix not keeping full alphabet in output
    for (a in setdiff(ALPHA_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM, vec))
    } 

    ##Selecting relevant columns
    Tidy_CM_Gaps <- select(Tidy_CM, ALPHA_21)
    Tidy_CM_NoGaps <- select(Tidy_CM, ALPHA_20)

    ##Entropy Calculation
    ent <- apply(Tidy_CM_Gaps, 1, entropy, unit="log2") %>% as_tibble()
    colnames(ent)<-paste0("Entropy_", HOG_ID)

    ##Entropy Calculation Ignoring Gaps
    entNG <- apply(Tidy_CM_NoGaps, 1, entropy, unit="log2") %>% as_tibble()
    colnames(entNG)<-paste0("EntropyNoGaps_", HOG_ID)

    ##Save fraction of invariant positions with/without gaps and number of highly variable positions (without gaps)
    frbl <- length(which(ent == 0))/nrow(ent)
    frblng <- length(which(entNG == 0))/nrow(entNG)
    nHVsites <- length(which(entNG > HSEC))                          ####KEY CUTOFF PARAMETER
    seq_num <- length(names(MinAli))
    stats <- c(HOG_ID , frbl, frblng, nrow(ent), nHVsites, seq_num)
    result_list <- list(ent, entNG, stats)
    return(result_list)
}

### Results presenting

cult_out <- cal_entropy(cult_msa, MinGapFraction, MinGapBlockWidth, Alph_21, Alph_20, HOG, hvSiteEntCutoff)
land_out <- cal_entropy(land_msa, MinGapFraction, MinGapBlockWidth, Alph_21, Alph_20, HOG, hvSiteEntCutoff)
wild_out <- cal_entropy(wild_msa, MinGapFraction, MinGapBlockWidth, Alph_21, Alph_20, HOG, hvSiteEntCutoff)
total_out <- cal_entropy(total_msa, MinGapFraction, MinGapBlockWidth, Alph_21, Alph_20, HOG, hvSiteEntCutoff)

#### Summary table

cult_stat <- c("Cultivar", cult_out[[3]])
land_stat <- c("Landrace", land_out[[3]])
wild_stat <- c("Wild", wild_out[[3]])
total_stat <- c("Total", total_out[[3]])

stat_list <- list(cult_stat, land_stat, wild_stat, total_stat)
stat_table <- as.data.frame(do.call(rbind, stat_list))
colnames(stat_table) <- c('Status', 'HOG', 'Invar_percent', 'Invar_percent_NG', 'MSA_length','Seq_Num')

write.csv(stat_table, file=paste0(args[2], "/", HOG,"_stat_table.csv"), na="", quote=TRUE)


#### Generate entropy table & Plot entropy value

plot_entropy <- function(DATA, STATUS){
    # set theme function 
    theme_custom <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
        base_rect_size = base_size/22) 
    {
        theme_bw(base_size = base_size, base_family = base_family, 
            base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black", size = rel(1)),
                  axis.ticks.length=unit(.15, "cm"),
                  legend.key = element_blank(), 
                  strip.background = element_rect(fill = "white", colour = "black", size = rel(2)), 
                  complete = TRUE,
    #               axis.ticks.x = element_blank(),
                  axis.text.x = element_text(size=16, angle=0, color="black", 
                                             vjust=1,hjust=1,
                                             margin = margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.text.y = element_text(size=16, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=16, angle=90, margin = margin(t = 0, r = 15, b = 0, l = 0)),
                  plot.title=element_text(family='', colour='red', size=20, hjust=0)
                 )
    }

    # tidy up data
    DATA <- as.data.frame(DATA)
    DATA$row_num <- seq.int(nrow(DATA))
    HOG <- names(DATA)[1]

    # set figure size
    options(repr.plot.width = 8, repr.plot.height = 4)
    p <- ggplot(DATA, aes(x = row_num, y= !!ensym(HOG))) +
                geom_line() + 
                theme_custom() +
                ggtitle(STATUS)
    result_list <- list(p, DATA)
    return(result_list)
}

# run plot_entropy
cult <- plot_entropy(cult_out[1], "Cultivar")
land <- plot_entropy(land_out[1], "Landrace")
wild <- plot_entropy(wild_out[1], "Wild")
cult_ng <- plot_entropy(cult_out[2], "Cultivar")
land_ng <- plot_entropy(land_out[2], "Landrace")
wild_ng <- plot_entropy(wild_out[2], "Wild")

# plot each subplot
cult_plt <- cult[[1]]
land_plt <- land[[1]]
wild_plt <- wild[[1]]
cult_ng_plt <- cult_ng[[1]]
land_ng_plt <- land_ng[[1]]
wild_ng_plt <- wild_ng[[1]]

# generate indexed_data
cult_data <- cult[[2]]
land_data <- land[[2]]
wild_data <- wild[[2]]
cult_ng_data <- cult_ng[[2]]
land_ng_data <- land_ng[[2]]
wild_ng_data <- wild_ng[[2]]

# generate integrated table
ent_table <- merge(cult_data, land_data, by="row_num") %>% 
             merge(wild_data, by="row_num") %>% 
             merge(cult_ng_data, by="row_num") %>% 
             merge(land_ng_data, by="row_num") %>% 
             merge(wild_ng_data, by="row_num")
colnames(ent_table) <- c("AA_order","ent_cultivar", "ent_landrace", "ent_wild", 
                         "ent_ng_cultivar", "ent_ng_landrace", "ent_ng_wild")

ent_table2 <- mutate(ent_table, ent_C2W=ent_cultivar-ent_wild) %>% 
                    mutate(ent_L2W=ent_landrace-ent_wild) %>% 
                    mutate(ent_C2L=ent_cultivar-ent_landrace) %>% 
                    mutate(ent_ng_C2W=ent_ng_cultivar-ent_ng_wild) %>% 
                    mutate(ent_ng_L2W=ent_ng_landrace-ent_ng_wild) %>% 
                    mutate(ent_ng_C2L=ent_ng_cultivar-ent_ng_landrace)

write.csv(ent_table2, file=paste0(args[2], "/", HOG,"_entropy_table.csv"), na="", quote=TRUE)

# plot the entropy difference
C2W <- plot_entropy(ent_table2$ent_C2W,"ent Cultvar - Wild")
L2W <- plot_entropy(ent_table2$ent_L2W,"ent Landrace - Wild")
C2L <- plot_entropy(ent_table2$ent_C2L,"ent Cultivar - Landrace")
ng_C2W <- plot_entropy(ent_table2$ent_ng_C2W,"entNG Cultivar - Wild")
ng_L2W <- plot_entropy(ent_table2$ent_ng_L2W,"entNG Landrace - Wild")
ng_C2L <- plot_entropy(ent_table2$ent_ng_C2L,"entNG Cultivar - Landrace")

C2W_plt <- C2W[[1]]
L2W_plt <- L2W[[1]]
C2L_plt <- C2L[[1]]
ng_C2W_plt <- ng_C2W[[1]]
ng_L2W_plt <- ng_L2W[[1]]
ng_C2L_plt <- ng_C2L[[1]]

# plot section
options(repr.plot.width = 16, repr.plot.height = 24)
fig <- plot_grid(cult_plt, cult_ng_plt, 
                 land_plt, land_ng_plt, 
                 wild_plt, wild_ng_plt,
                 C2W_plt, ng_C2W_plt,
                 L2W_plt, ng_L2W_plt,
                 C2L_plt, ng_C2L_plt,
                 ncol=2,
                 scale=0.9)

ggsave2(file=paste0(args[2], "/", HOG,".png"), fig, height=24, width=16, dpi=300)
