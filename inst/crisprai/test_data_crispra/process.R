# TSS
tss <- read.table("crispra_tss_example.txt", head=TRUE)
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
grna <- read.table("crispra_sgrna_example.txt", head=TRUE)
grna <- grna[grna$tss_id %in% tss$tss_id,]



tssExampleCrispra  <- tss
sgrnaExampleCrispra <- grna
save(tssExampleCrispra,
     file="../../../data/tssExampleCrispra.rda")
save(sgrnaExampleCrispra,
     file="../../../data/sgrnaExampleCrispra.rda")
