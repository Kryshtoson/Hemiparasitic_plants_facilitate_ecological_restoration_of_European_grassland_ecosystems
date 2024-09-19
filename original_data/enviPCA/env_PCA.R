library(vegan)
library(psych)

# envi <- read.delim("clipboard", row.names = 1)
envi <- read.csv("enviPCA/data/sites_with_soil_JT.csv", row.names = 1,
               skip=1)
envi <- envi[,1:11]

pdf("envi.pairs.pdf", 8,8)
pairs.panels(envi)
dev.off()
names(envi)

envi[,c(6,9)] <- log(envi[,c(6,9)])

pca.1<-rda(envi, scale=T)
pca.1
screeplot(pca.1, bstick=T)



aa<-ordiplot(pca.1, display=c("sp", "si"), type="n", scaling=3)
text(aa, what="sp")
text(aa, what="si", col=2)


aa<-ordiplot(pca.1, display=c("sp", "si"), type="n", scaling=3, choices=c(2, 3))
text(aa, what="sp")
text(aa, what="si", col=2)
pcs<-data.frame(scores(pca.1, choices=1:3, display="si"))
pcs$site=rownames(pcs)
