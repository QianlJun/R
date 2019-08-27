library(biomaRt)
listMarts()
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
class(mouse)
genes = markers
genes = getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
               values = genes, mart = mouse,
               attributesL = c("hgnc_symbol","chromosome_name","start_position"),
               martL = human,
               uniqueRows = T)
genes
