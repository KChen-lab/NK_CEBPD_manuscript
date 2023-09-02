
##########   Fig 4C Motif enrichment for cluster 6 specific accessibility regions ########


# compare of TFs in cluster 2 and 5 
DefaultAssay(part_rmIL12) <- "peaks"
da_peaks <- FindMarkers(
  object = part_rmIL12,
  min.pct = 0.2,
  only.pos = FALSE,
  ident.1=c("5"),
  ident.2=c("2"),
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
   latent.vars = "nCount_peaks",
   logfc.threshold=0.05
  )
up_peaks <- rownames(da_peaks)[da_peaks$avg_diff<0]

clst5_motifs <- FindMotifs(object = part_rmIL12,features = up_peaks)


clst5_motifs$pvalue[clst5_motifs$pvalue<10^(-200)]<-10^(-200)

clst5_motifs$fold.enrichment[grepl("CEB", clst5_motifs$motif.name)] <- 1 + clst5_motifs$fold.enrichment[grepl("CEB", clst5_motifs$motif.name)]
plt_dt <- data.frame("fold_enrichment"=clst5_motifs$fold.enrichment, 
  "pvalue"=(0-log10(clst5_motifs$pvalue)),
  "name"=clst5_motifs$motif.name)

#plt_dt$name[plt_dt$fold_enrichment<=3]<-""
plt_dt$name[!plt_dt$name%in%c("EOMES","NFKB2","IRF1","STAT1","CEBPD","ATF4","CEBPB","JUN","RUNX2","KLF13","TBX21","KLF2",
  "KLF6","FLI1","REST","E2F3","ZEB1","STAT2","FOXO3","MAX","SIN3A","MXI1","RELB","STAT3","ETS1","IRF9")] <-""

p <- ggplot(plt_dt, aes(fold_enrichment, pvalue)) +
  geom_point(color = "gray") + theme_bw() +  
  geom_point(data=plt_dt[!plt_dt$name=="",], aes(fold_enrichment, pvalue), color="#FF7F50") + 
  geom_text_repel(aes(label = plt_dt$name),color="#FF7F50",size = 3.5, max.overlaps=Inf) + xlab("Fold (Enrichment)") + 
  ylab("-log10(pval)") + ylim(0,250)

tiff("./atac_umap_clst5_TF.tiff", width=6, height=6, res =300, units = "in", compression = "lzw")
print(p)
dev.off()

