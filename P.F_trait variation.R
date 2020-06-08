




#   EXTRACT DATA FROM GITHUB  

library(curl)
All_traits <- read.csv( curl("https://raw.githubusercontent.com/dblasini/FunctTraits/master/P.F_traits.csv") )
attach (All_traits )

All_traits 




#########################               LEAF SPECTRUM   


# ONLY WITH ECOTYPE AND ELEVATION 
Leaf.traits<-All_traits[, c("Ecotype","Elevation","Ds", "Ss","gsmax",
                            "SLA","LS","Ap","PF","PLuA","PHMD","PLuF","PVD","Kh","Kl")]

# ONLY WITH ECOTYPE  AND GENOTYPE
Leaf.traits.2<-All_traits[, c("Genotype","Elevation","Ds", "Ss","gsmax",
                             "SLA","LS","Ap","PF","PLuA","PHMD","PLuF","PVD","Kh","Kl")]

# HERE, WE RUN PCR 

# BECAUSE THE EXTREME HIGH CORRELATION BETWEEN gsmax AND Ds,
# I DECIDED TO ELIMINATE THIS TRAIT FROM THE PCA ANALYSIS 

traits.L <- prcomp(Leaf.traits[c(1:48),c(-1:-2, -5)], center=TRUE, scale=TRUE)
traits.L

summary(traits.L)

library(vegan)
screeplot(traits.L,bstick=TRUE)


# LINEAR REGRESSION BETWEEN LEAF PC1 AND THE 8 PROVENANCE ELEVATIONS 
modelpc1<-lm(Leaf.traits[,-1]$Elevation ~ traits.L$x[,1], data = Leaf.traits)
summary(modelpc1)    #  p-value: 1.432e-08, Adjusted R-squared:  0.4956 

#   FIRST AXIS PLOT
par(mfrow=c(1,1))
plot(Leaf.traits[,-1]$Elevation, traits.L$x[,1], yaxt='n', ann=FALSE)
legend("topleft", cex=1,legend=c(expression(R^2 == 0.50, "p < 0.001")), pch=c())
title(main="", xlab="Elevation (m)", ylab="Leaf Trait PC1",cex.lab=1.5)




# LINEAR REGRESSION BETWEEN LEAF PC2 AND THE 8 PROVENANCE ELEVATIONS 
modelpc2<-lm(Leaf.traits[,-1]$Elevation ~ traits.L$x[,2], data = Leaf.traits)
summary(modelpc2)    #  p-value: 0.004633, Adjusted R-squared:  0.1433 
#   SECOND AXIS PLOT
par(mfrow=c(1,1))
plot(Leaf.traits[,-1]$Elevation, traits.L$x[,2], yaxt='n', ann=FALSE)
legend("topright", cex=1,legend=c(expression(R^2 == 0.16, "p < 0.01")), pch=c())
title(main="", xlab="Elevation (m)", ylab="Leaf Trait PC2",cex.lab=1.5)




###      PCA visualization

library(tidyverse)
library(magrittr)

#   Principal Component Analysis
library("FactoMineR")
library("factoextra")


L.pca <- PCA(Leaf.traits [c(1:48),c(-1, -5)], graph = FALSE)
eig.val <- get_eigenvalue(L.pca )
eig.val


fviz_eig(L.pca, addlabels = TRUE, ylim = c(0, 50))

# Contributions of variables to PC1
fviz_contrib(L.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(L.pca, choice = "var", axes = 2, top = 10)

# Description of dimension 1
LT.PC1 <- res.desc$Dim.1
LT.PC1

# Description of dimension 2
LT.PC2 <- res.desc$Dim.2
LT.PC2



Leaf.Traits.Ecotype<- fviz_pca_biplot(L.pca, 
                       geom.ind = "point",
                       col.ind = Leaf.traits$Ecotype,
                       pointsize = 4,labelsize = 6,
                       palette = "npg", 
                       addEllipses = TRUE, label = "var", # "none" make no labels 
                       col.var = "black", repel = TRUE,
                       legend.title = "Ecotype") 
Leaf.Traits.Ecotype + theme(plot.title = element_text(color="black", size=18, face="bold"),
             axis.title.x = element_text(color="black", size=15, face="bold"),
             axis.title.y = element_text(color="black", size=15, face="bold"))








####################     ARCHITECTURE CORNER'S RULE SPECTRUM 



# ONLY WITH ECOTYPE AND ELEVATION 
Architecture.traits<-All_traits[c(1:48), c(1,4,20:24)]

# ONLY WITH ECOTYPE AND GENOTYPE
Architecture.traits.2<-All_traits[c(1:48), c(1,3,20:24)]


# HERE, WE RUN PCR 
traits.A <- prcomp(Architecture.traits[c(1:48),c(-1:-2)], center=TRUE, scale=TRUE)
traits.A

summary(traits.A)

library(vegan)
screeplot(traits.A,bstick=TRUE)


# LINEAR REGRESSION BETWEEN LEAF PC1 AND THE 8 PROVENANCE ELEVATIONS 
modelpc1A<-lm(Architecture.traits[,-1]$Elevation ~ traits.A$x[,1], data = Architecture.traits)
summary(modelpc1A)    #  p-value: 3.141e-09, Adjusted R-squared:  0.5273 



###      PCA visualization

library(tidyverse)
library(magrittr)

#   Principal Component Analysis
library("FactoMineR")
library("factoextra")


A.pca <- PCA(Architecture.traits [c(1:48),c(-1:-2)], graph = FALSE)
eig.val <- get_eigenvalue(A.pca )
eig.val


fviz_eig(A.pca, addlabels = TRUE, ylim = c(0, 50))

# Contributions of variables to PC1
fviz_contrib(A.pca, choice = "var", axes = 1, top = 10)



Architecture.Traits.Ecotype<- fviz_pca_biplot(A.pca, 
                                      geom.ind = "point",
                                      col.ind = Architecture.traits$Ecotype,
                                      pointsize = 4,labelsize = 6,
                                      palette = "npg", 
                                      addEllipses = TRUE, label = "var", # "none" make no labels 
                                      col.var = "black", repel = TRUE,
                                      legend.title = "Ecotype") 
Architecture.Traits.Ecotype + theme(plot.title = element_text(color="black", size=18, face="bold"),
                            axis.title.x = element_text(color="black", size=15, face="bold"),
                            axis.title.y = element_text(color="black", size=15, face="bold"))









