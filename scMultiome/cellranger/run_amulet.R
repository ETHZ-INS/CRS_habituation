library(scDblFinder)

lf <- paste0(list.files(pattern="sample_"),"/outs")
for(f in lf){
message(f)
s <- amulet(paste0(f,"/atac_fragments.tsv.gz"), barcodes=paste0(f,"/cellbender_filtered_cell_barcodes.csv"))
saveRDS(s, file=paste0(f,"/atac_amulet_results.rds"))	
}

