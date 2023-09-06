### create gene count file for each sample ###
library(Rsubread)
BAMs <- list.files(pattern="*.bam")
# inbuilt has been changed to RefSeqGTF UCSC 2020
for (f in 1:length(BAMs)){ 
  fc_mx <- featureCounts(BAMs[f], annot.inbuilt = "hg38", 
                        isPairedEnd = TRUE, countChimericFragments = FALSE, 
                        requireBothEndsMapped = TRUE, countMultiMappingReads = FALSE)
  write.csv(fc_mx$counts, paste(strsplit(BAMs[f], "[_]")[[1]][1], ".csv", sep=""))
}


### normalize gene count with rlog method ###
library(DESeq2)
fcs <- list.files(pattern="*.csv")
fcm <- read.table(fcs[1], sep=",", header = TRUE, stringsAsFactors = FALSE)
names(fcm)[1] <- "geneid"
names(fcm)[2] <- strsplit(fcs[1], "[.]")[[1]][1]
for (f in 2:length(fcs)){
  fcm2 <- read.table(fcs[f], sep=",", header = TRUE, stringsAsFactors = FALSE)
  fcm <- cbind(fcm, fcm2[, 2])
  names(fcm)[f+1] <- strsplit(fcs[f], "[.]")[[1]][1]
}
row.names(fcm) <- fcm$geneid
fcm[1] <- NULL
coldata <- data.frame(colnames(fcm))
# save gene count matrix
write.csv(fcm, "results/count_matrix.csv")

# transform with rlog method, save log scale relative expression matrix
dds <- DESeqDataSetFromMatrix(countData = fcm, colData = coldata, design = ~ 1)
logm <- rlog(dds, blind=TRUE)
write.csv(assay(logm), "results/rlogm.csv")


### compare two groups ###
# AMnp day3 vs day14
coldata <- data.frame(colnames(fcm[28:36]))
dds <- DESeqDataSetFromMatrix(countData = fcm[28:36], colData = coldata, design = ~ 1)
dds$day <- factor(rep(rep(c("3", "7", "14")), 3))
design(dds) <- ~ day
DEm <- DESeq(dds, test = "Wald", fitType = "parametric")
result_Amnp_3_14 <- results(DEm, format = "DataFrame", contrast = c("day", "3", "14"))
write.csv(result_Amnp_3_14, "results/Amnp_day3_day14.csv")


### volcano plot for comparing two groups ###
library(EnhancedVolcano)
library(tidyverse)
library(cowplot)
library(ggsci)
library(ggrepel)
library(ggplot2)
library(ggthemes)

my_palette <- c( "#4DBBD5FF", "#999999", "#E64B35FF")

rbps <- c('CELF5', 'CPEB3', 'ELAVL4', 'RBFOX3',
          'ELAVL2', 'LIN28B', 'KHDRBS2')
gaba <- c('DLX1', 'DLX2', 'DLX5','DLX6','LHX6','GBX1', 'GBX2', 'OTX2',
          'LAMP5','VIP','GAD1','GAD2','SLC6A1','SLC17A8','SLC32A1')
fib <- c('DCN','FN1','FBLN1','FBLN2','DMKN','HSPB1','PTX3','TGFBI','COL1A2','COL13A1',
         'S100A6','S100A4','LGALS1','LOXL1','KRT5')

Amnp_D3_14 <- read_csv("Amnp_day3_day14.csv")
names(Amnp_D3_14)[2] <- "geneid"
# remove genes with NA value
Amnp_D3_14_de <- na.omit(Amnp_D3_14_de)

# select differentially expressed genes & annotation
group_by(Amnp_D3_14_de, direction) %>%
  summarise(n())
Amnp_D3_14_de <- select(Amnp_D3_14, geneid, log2FoldChange, padj) %>%
  mutate(FC = 2**log2FoldChange, 
         direction = if_else(padj > 0.05 | abs(log2FoldChange) < 1, 'ns',
                             if_else(log2FoldChange >= 1, "up", "down")))%>%
  mutate(selected = if_else(geneid %in% gaba, 
                            'Yes', direction))
# make plot
ggplot(data = Amnp_D3_14_de, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 3, aes(color = direction)) +
  geom_point(data = filter(Amnp_D3_14_de, selected == 'Yes'),
             shape = 21, 
             size = 3, 
             color = 'black',
             stroke = 0.8)+
  geom_text_repel(data = filter(Amnp_D3_14_de, selected == 'Yes'),
                  size = 5, 
                  nudge_x = .10,
                  box.padding = 0.15,
                  aes(label = geneid),
                  max.overlaps = 40) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotdash',
             color = 'grey30') +
  geom_vline(xintercept = c(-1, 1),
             linetype = 'dotdash',
             color = 'grey30') +
  scale_color_manual(values = my_palette)+
  ylim(0, 400) +
  xlim(-12, 12) +
  labs(x = 'Log2(Fold Change)',
       y = '-log10(p.adjust)',
       title = "AMnp (day14_day3)",
       size = "log2 fold change",
       color = "") +
  theme_few()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = c(0.1, 0.88),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
)


### gene ontology plot ###
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(DOSE)
library(ggsci)
library(cowplot)
library(tidyverse)
library(factoextra)
library(AnnotationDbi)
library(clusterProfiler) # version 4.5.2
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(ggpubr)

# enrichGO_Amnp_DEG_up(493 genes)
Amnp_DEG_up <- read_csv("Amnp_DEG_up_gene.ID.csv")

# BP annotation
Amnp_DEG_up_enrich <- enrichGO(Amnp_DEG_up$ENTREZID, 
                               OrgDb = "org.Hs.eg.db",
                               keyType = "ENTREZID", 
                               ont = "ALL",
                               readable = T,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 5000)

Amnp_DEG_up_enrich <- as.data.frame(Amnp_DEG_up_enrich)

write.csv(Amnp_DEG_up_enrich, "Amnp_DEG_up_enrich.csv", row.names = F)

# enrichGO_Amnp_DEG_down(104 genes)
Amnp_DEG_down <- read_csv("RNA-seq/venn/Amnp_DEG_down_gene.ID.csv")

# BP annotation
Amnp_DEG_down_enrich <- enrichGO(Amnp_DEG_down$ENTREZID, 
                                 OrgDb = "org.Hs.eg.db",
                                 keyType = "ENTREZID", 
                                 ont = "ALL",
                                 readable = T,
                                 minGSSize = 1,
                                 maxGSSize = 5000)

Amnp_DEG_down_enrich <- as.data.frame(Amnp_DEG_down_enrich)

write.csv(Amnp_DEG_down_enrich, file = "Amnp_DEG_down_enrich.csv")

# bar plots
p1 <- arrange(Amnp_DEG_up[1:10, ])
p2 <- p1 %>% mutate(log10 = -log10(p.adjust))
p2 <- arrange(p2, log10)
p2$Description <- factor(p2$Description, levels=p2$Description)

ggplot(p2, aes(x=-log10(p.adjust), y=Description))+
    geom_bar(stat="identity", width=0.8, aes(color = log10(p.adjust)), fill = "#E64B35")+
    xlab("-log10(p.adjust)")+
    ylab("")+
    ggtitle("Biological Process (AMnp_up)")+
    scale_x_continuous(limits = c(0, 40))+
    theme_bw()+
    theme(axis.title.x = element_text(size = 22, color="black", family="Arial"),
          axis.title.y = element_blank())+
    theme(axis.text.x = element_text(size = 22, color="black", family="Arial"),
        axis.text.y = element_text(size = 22, color="black", family="Arial"),
        legend.position="none",
        plot.title = element_text(size = 22)
)



p1 <- arrange(Amnp_DEG_down_2[1:10, ])
p2 <- p1 %>% mutate(log10 = -log10(p.adjust))
p2 <- arrange(p2, desc(log10(p.adjust)))
p2$Description <- factor(p2$Description, levels=p2$Description)

ggplot(p2, aes(x=-log10(p.adjust), y=Description))+
  geom_bar(stat="identity", width=0.8, aes(color = log10(p.adjust)), fill = "#4DBBD5")+
  xlab("-log10(p.adjust)")+
  ylab("")+
  ggtitle("Biological Process (AMnp_down)")+
  scale_x_continuous(limits = c(0, 40))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 21, color="black", family="Arial"),
        axis.title.y = element_blank())+
  theme(axis.text.x = element_text(size = 21, color="black", family="Arial"),
        axis.text.y = element_text(size = 21, color="black", family="Arial"),
        legend.position="none",
        plot.title = element_text(size = 21)
)