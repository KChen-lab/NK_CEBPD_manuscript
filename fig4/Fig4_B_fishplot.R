

############  Fig 4B Dynamical cluster change of for scATAC-seq: IL12, IL15 on 3 timepoints ########
library(ggplot2)
library(reshape2)
library(ggsci)
library(timescape)
options(browser="/usr/bin/firefox")
library(reshape2)

sta <- table(part_rmIL12@meta.data[,c("seurat_clusters","product")])
for(i in seq(1,nrow(sta),1)){
  sta[i,] <- sta[i,]/sum(sta[i,])
}

sta<-round(sta,digits=3)
sta <- sta*100
sta <- melt(sta)
sta <- as.data.frame(sta)

colnames(sta)<-c("Cluster","Product","Percentage")
sta$label <- sta$Percentage
sta$label[sta$label<2] <- ""
sta$Cluster <- as.factor(sta$Cluster)
sta$Product <- factor(sta$Product, levels=c("IL21_GBM_9","IL15_GBM_9","IL21_GBM_3","IL15_GBM_3","IL21_baseline","IL15_baseline"))


prev <- melt(atac_sta)
colnames(prev) <- c("clone_id","timepoint","clonal_prev")
prev$clonal_prev[is.na(prev$clonal_prev)]<-0
edges <- data.frame("source" = rep("N0", nrow(sta)),
                    "target" =  as.character(unique(prev$clone_id))[seq(1, nrow(sta))])

tp<-prev[seq(1, nrow(sta),1),]
tp$clone_id<-"N0"
tp$clonal_prev <- 0

prev<-rbind(prev,tp)

#mycolor_code <- col2hex(nk_color)

#clone_color <- data.frame(clone_id = as.character(unique(prev$clone_id)), colour =c( mycolor_code,"#808080"))
clone_color <- data.frame(clone_id = c(seq(0,nrow(sta),1)),
                          colour =c( vega_20_scanpy[seq(1, nrow(sta),1)],"#808080"))
nk_plt<-timescape(clonal_prev =prev[grepl("IL15", prev[,2]), ], tree_edges = edges, clone_colours =clone_color , height=400, alpha=15)
htmlwidgets::saveWidget(nk_plt,file.path(paste0("~/project/nk_multiOme/atac.NK.IL15.html")))
nk_plt<-timescape(clonal_prev =prev[grepl("IL21", prev[,2]), ], tree_edges = edges, clone_colours =clone_color , height=400, alpha=15)
htmlwidgets::saveWidget(nk_plt,file.path(paste0("~/project/nk_multiOme/atac.NK.IL21.html")))

