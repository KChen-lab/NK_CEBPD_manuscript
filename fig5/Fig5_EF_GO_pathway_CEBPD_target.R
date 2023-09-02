

################## Fig 5 E-F GO pathway for CEBPB and CEBPD target genes  ###################

library(reshape2)
library(ggplot2)


dt <- read.csv(file="./../result_Jan1_2022/atac/CEBP_TFs/CEBPD_pathway.csv")

dt$Count[7]<-7
dt$adj_pval <- (0-log10(dt$P.value))
dt$Term <- factor(dt$Term, levels=dt$Term)
dt$Count <- as.numeric(as.character(dt$Count))

dt <- melt(dt[,c('Term','adj_pval','Count')],id.vars = 1)
dt$value[dt$variable=="Count"] <- (0-dt$value[dt$variable=="Count"])
p <- ggplot(dt,aes(x = Term,y = value)) + 
    geom_bar(aes(fill = variable),stat = "identity",position = "identity")  + 
    labs(x=NULL, y=expression(-log[10] * ' p-value'))+
    coord_flip()+ theme_classic() + 
    scale_fill_manual(values = c("burlywood1", "steelblue1"))


pdf(file="./CEBPD_GO.pdf",height=4,width=6)
print(p)
dev.off()




dt <- read.csv(file="./../result_Jan1_2022/atac/CEBP_TFs/CEBPB_HallmarkPathway.csv")

dt$Count[7]<-7
dt$adj_pval <- (0-log10(dt$P.value))
dt$Term <- factor(dt$Term, levels=dt$Term)
dt$Count <- as.numeric(as.character(dt$Count))

dt <- melt(dt[,c('Term','adj_pval','Count')],id.vars = 1)
dt$value[dt$variable=="Count"] <- (0-dt$value[dt$variable=="Count"])
p <- ggplot(dt,aes(x = Term,y = value)) + 
    geom_bar(aes(fill = variable),stat = "identity",position = "identity")  + 
    labs(x=NULL, y=expression(-log[10] * ' p-value'))+
    coord_flip()+ theme_classic() + 
    scale_fill_manual(values = c("burlywood1", "steelblue1"))

pdf(file="./CEBPB_GO.pdf",height=2.5,width=6)
print(p)
dev.off()

