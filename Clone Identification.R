rm(list = ls())
setwd("")

# generate color 
genCol <- function(x, trans){
 
 return(rgb(red = col2rgb(x, alpha = FALSE)[1],
            green = col2rgb(x, alpha = FALSE)[2],
            blue = col2rgb(x, alpha = FALSE)[3],
            max = 255, alpha = trans, names = NULL))
 
}


library("WGCNA") 
library("sparcl")
library("data.table")

# reading IBS matrix based on SNPs with allele frequency >= 0.05:

maAll = as.matrix(read.table("results.ibsMat"))

maBlue = as.matrix(read.table("results1.ibsMat"))

maBrown = as.matrix(read.table("results2.ibsMat"))

maRed = as.matrix(read.table("results3.ibsMat"))

bamBlue <- readLines("results1.ibsMat")


bamBrown <- readLines("results2.ibsMat")

bamRed <- readLines("results3.ibsMat")

fread("results.ibsMat")
resAll <- fread("results.ibsMat")

fread("results.ibsMat")
res1 <- fread("results1.ibsMat")

fread("results.ibsMat")
res2 <- fread("results2.ibsMat")

fread("results.ibsMat")
res3 <- fread("results3.ibsMat")

#### All dendogram



# plotting hierarchical clustering tree

hcAll=hclust(as.dist(maAll),"ave")
plot(hcAll,cex=0.6)

# sorting samples into clonal groups; 
# singletons go into the same "color" group (for plotting) 
# but different "cn" groups (for GLM modeling later)
cc=cutree(hcAll, h = 0.003322)

# meta$cn=cc

tc=table(cc)

singletons=as.numeric(names(tc)[tc==1])

cc[cc %in% singletons]= singletons[1]

clones=labels2colors(as.numeric(as.factor(as.numeric(cc))))

clones[cc == 1] <- "black"

pdf("DendroAll.pdf", height = 8, width = 12)

ColorDendrogram(hclust(as.dist(maAll),"ave"), y = clones, 
                branchlength=0.035, xlab = "All populations", main = "", ylab = "Genetic Distance (IBS)",
                labels = rownames(resAll))


dev.off()


#### Population 1 dendogram


# plotting hierarchical clustering tree
hc=hclust(as.dist(maBlue),"ave")
plot(hc,cex=0.6)

# sorting samples into clonal groups; 
# singletons go into the same "color" group (for plotting) 
# but different "cn" groups (for GLM modeling later)
cc=cutree(hc, h = 0.0032836)

# meta$cn=cc

tc=table(cc)

singletons=as.numeric(names(tc)[tc==1])

cc[cc %in% singletons]= singletons[1]

clones=labels2colors(as.numeric(as.factor(as.numeric(cc))))

clones[cc == 1] <- "black"

pdf("Dendro1.pdf", height = 8, width = 12)

ColorDendrogram(hclust(as.dist(maBlue),"ave"), y = clones, 
                branchlength=0.035, xlab = "Population 1", main = "", ylab = "Genetic Distance (IBS)",
                labels = rownames(res1))

#### Population 2 dendogram


# plotting hierarchical clustering tree
hc=hclust(as.dist(maBrown),"ave")
plot(hc,cex=0.6)

# sorting samples into clonal groups; 
# singletons go into the same "color" group (for plotting) 
# but different "cn" groups (for GLM modeling later)
cc=cutree(hc, h = 0.0033258)

# meta$cn=cc

tc=table(cc)

singletons=as.numeric(names(tc)[tc==1])

cc[cc %in% singletons]= singletons[1]

clones=labels2colors(as.numeric(as.factor(as.numeric(cc))))

clones[cc == 1] <- "black"

pdf("Dendro1.pdf", height = 8, width = 12)

ColorDendrogram(hclust(as.dist(maBrown),"ave"), y = clones, 
                branchlength=0.035, xlab = "Population 2", main = "", ylab = "Genetic Distance (IBS)",
                labels = rownames(res2))


#### Population 3 dendogram


# plotting hierarchical clustering tree
hc=hclust(as.dist(maRed),"ave")
plot(hc,cex=0.6)

# sorting samples into clonal groups; 
# singletons go into the same "color" group (for plotting) 
# but different "cn" groups (for GLM modeling later)
cc=cutree(hc, h = 0.0033258)

# meta$cn=cc

tc=table(cc)

singletons=as.numeric(names(tc)[tc==1])

cc[cc %in% singletons]= singletons[1]

clones=labels2colors(as.numeric(as.factor(as.numeric(cc))))

clones[cc == 1] <- "black"

pdf("Dendro1.pdf", height = 8, width = 12)

ColorDendrogram(hclust(as.dist(maRed),"ave"), y = clones, 
                branchlength=0.035, xlab = "Population 3", main = "", ylab = "Genetic Distance (IBS)",
                labels = rownames(res3))
