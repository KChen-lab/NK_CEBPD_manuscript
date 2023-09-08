#################### Fig 5 C-D VlnPlot for CEBPB and CEBPD difference between IL21 and IL15 ###################

figs <- readRDS(file="CEBPD.target.RDS")
p <- figs
p$data$id <- as.character(p$data$id)

mycompare <- list(c("IL15","IL21"))
p$data$avg.exp.scaled[p$data$avg.exp.scaled>1]<-1
p$data$avg.exp.scaled[p$data$avg.exp.scaled<(-1)]<-(01)
plt_dt <- p$data[p$data$product%in%c("IL15","IL21"),]
plt_dt$product <- factor(plt_dt$product,levels=c("IL21","IL15"))
p1 <- ggplot(plt_dt,aes(x=product, y=avg.exp.scaled)) + 
  geom_boxplot(aes(fill=product))  + facet_wrap(~time)  + 
 theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))  + 
 ylab("Scaled gene-level chromatin accessibility ") + stat_compare_means(comparisons = mycompare)

p2 <- ggplot(plt_dt,aes(x=product, y=avg.exp.scaled)) + 
  geom_violin(aes(fill=product))  + facet_wrap(~time)  + 
 theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))  + 
 ylab("Scaled gene-level chromatin accessibility ") + stat_compare_means(comparisons = mycompare)

pdf("./CEBPD.target.pdf", width=6, height =4)
print(p1)
print(p2)
dev.off()
