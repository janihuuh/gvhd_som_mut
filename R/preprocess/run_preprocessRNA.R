

message("Create individual seurat-objects ...")

t2013 <- Read10X(data.dir = "data/scRNAseq/2013/")  %>% CreateSeuratObject(project = "2013", min.cells = 3, min.features = 200)
t2017 <- Read10X(data.dir = "data/scRNAseq/2017/")  %>% CreateSeuratObject(project = "2017", min.cells = 3, min.features = 200)

## Merge time points 
gvhd   <- merge(t2013, t2017, add.cell.ids = c("2013", "2017"))

## Basic QC
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

gvhd  <- PercentageFeatureSet(gvhd, pattern = "^MT-", col.name = "percent.mt")
gvhd  <- PercentageFeatureSet(gvhd, pattern = "^RP", col.name = "percent.ribo")
gvhd  <- PercentageFeatureSet(gvhd, features = cycle.genes, col.name = "percent.cycle")

gvhd@meta.data$barcode   <- colnames(gvhd)
gvhd@meta.data$timepoint <- substr(colnames(gvhd), 1, 4)
