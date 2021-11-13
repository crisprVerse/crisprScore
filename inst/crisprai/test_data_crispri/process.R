# TSS
tss <- read.table("crispri_tss_example.txt", head=TRUE)
genesToKeep <- c("A1BG", "SMARCA2", "KRAS")
tss <- tss[tss$gene_symbol %in% genesToKeep,]
tss$tss_id <- paste0(tss$gene_symbol, "_", tss$promoter)
rownames(tss) <- NULL
cols <- c("tss_id",
          "gene_symbol",
          "promoter",
          "transcripts",
          "position",
          "strand",
          "chr")
tss <- tss[,cols]

# GRNA:
grna <- read.table("crispri_sgrna_example.txt", head=TRUE)
grna <- grna[grna$tss_id %in% tss$tss_id,]



tssExampleCrispri  <- tss
sgrnaExampleCrispri <- grna
save(tssExampleCrispri,
     file="../../../data/tssExampleCrispri.rda")
save(sgrnaExampleCrispri,
     file="../../../data/sgrnaExampleCrispri.rda")