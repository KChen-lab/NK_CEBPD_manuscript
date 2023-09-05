
###########   Fig 4D Motif-based TF activity level difference between IL21 and Il15 ########

#day3$pvalue<-log10(day3$pvalue)
tf_change <- readRDS(file="./tf_day3_9_il15_21.RDS")
day3<-tf_change$day3


day9 <- tf_change$day9


day3$pvalue[day3$pvalue<10^(-20)]<-10^(-20)

p1<-  EnhancedVolcano(day3,
    lab = day3$motif.name,
    x = 'fold.enrichment',
    subtitle="", widthConnectors=1,
    y = 'pvalue',
    title = 'IL15_day3 versus IL21_day3',
    pCutoff = 10e-5,
    FCcutoff = 2,
    pointSize = 1.0,
    labSize = 3.0,
    col=c('gray', 'gray', 'gray', 'red'),
    colAlpha = 1)
day9$pvalue[day9$pvalue<10^(-20)]<-10^(-20)

p2 <- EnhancedVolcano(day9,
    lab = day9$motif.name,
    x = 'fold.enrichment',
    y = 'pvalue',
    title = 'IL15_day9 versus IL21_day9', subtitle="", widthConnectors=0.5,
    pCutoff = 10e-5,
    FCcutoff = 2,
    pointSize = 1.0,
    labSize = 3.0,
    col=c('gray', 'gray', 'gray', 'red'),
    colAlpha = 1)



tiff("TF_il15_21_day3.tiff", width=7, height =7, units = 'in', res=300, compression="lzw")
#DoHeatmap(obj$obj, features=sigFeature$gene,slot="data")
print(p1)
dev.off()


tiff("TF_il15_21_day9.tiff", width=7, height =7, units = 'in', res=300, compression="lzw")
#DoHeatmap(obj$obj, features=sigFeature$gene,slot="data")
print(p2)
dev.off()
