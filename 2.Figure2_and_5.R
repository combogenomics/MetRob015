require(pheatmap)
#library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggplot2)
library(reshape2)
library(heatmap3)
library(forcats)
library(grid)
library(gridExtra)
library('corrplot')
library(RColorBrewer)
library(psych)
library(vegan)
library(geosphere)
library(metan)
library(MatrixCorrelation)
library(cultevo)
library(ade4)
library(LambertW)
library(matrixTests)
library(jaccard)
library(qvalue)


wd <- '/Users/marcofondi/Dropbox/PhTAC125/AnTReN/metabolomics/'
setwd(wd)

#import the data
zerodata1 <- read.table("Intra_0.txt", sep = '\t', header = T, dec = ',')

#average every two lines
zerodata2 <- do.call(rbind,
        lapply(seq(1, nrow(zerodata1), 2), function(i){
          x <- zerodata1[ i:(i + 1), , drop = FALSE]
          res <- rbind( colSums(x)/2)
          rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
          res
        }))

#compute sd every two lines
zerostd_df <- do.call(rbind,
                 lapply(seq(1, nrow(zerodata1), 2), function(i){
                   x <- zerodata1[ i:(i + 1), , drop = FALSE]
                   res <- rbind( apply(x, 2, sd, na.rm = TRUE))
                   rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                   res
                 }))
#convert into dataframe
zerodata3 <- data.frame(zerodata2)
rownames(zerodata3)<-zerodata3$OD
rownames(zerodata3) <- as.numeric(c(1,2,3,4,5))

#normalize by OD
#data4 <- data3[,2:ncol(data3)]/as.numeric(rownames(data3))
zerodata4 <- zerodata3
zerodata4$OD <- zerodata3$OD
#order by OD
zerodata5 <- zerodata4
zerodata6 <- as.matrix(zerodata5[,2:(ncol(zerodata5))])


zerodata5_m <- melt(t(zerodata5[,2:(ncol(zerodata5))]), id.vars = as.numeric(zerodata5$time))
zerosd_m <- melt(t(zerostd_df[,2:(ncol(zerostd_df))]), id.vars = as.numeric(zerostd_df$time))
zerodata5_m$sd <- zerosd_m$value


pcold <- ggplot(data=zerodata5_m, aes(x=Var2, y=value, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point(aes(color=Var1))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=Var1),  width=.1, position=position_dodge(0.05))+
  facet_grid(Var1 ~ ., scales="free")+
    theme_classic()
pcold


#import 15 degrees data

#import the data
fifteendata1 <- read.table("Intra_15.txt", sep = '\t', header = T, dec = ',')

#average every two lines
fifteendata2 <- do.call(rbind,
                     lapply(seq(1, nrow(fifteendata1), 2), function(i){
                       x <- fifteendata1[ i:(i + 1), , drop = FALSE]
                       res <- rbind( colSums(x)/2)
                       rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                       res
                     }))

#compute sd every two lines
fifteenstd_df <- do.call(rbind,
                      lapply(seq(1, nrow(fifteendata1), 2), function(i){
                        x <- fifteendata1[ i:(i + 1), , drop = FALSE]
                        res <- rbind( apply(x, 2, sd, na.rm = TRUE))
                        rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                        res
                      }))
#convert into dataframe
fifteendata3 <- data.frame(fifteendata2)
rownames(fifteendata3)<-fifteendata3$OD
rownames(fifteendata3) <- as.numeric(c(1,2,3,4,5))

#normalize by OD
#data4 <- data3[,2:ncol(data3)]/as.numeric(rownames(data3))
fifteendata4 <- fifteendata3
fifteendata4$OD <- fifteendata3$OD
#order by OD
fifteendata5 <- fifteendata4
fifteendata6 <- as.matrix(fifteendata5[,2:(ncol(fifteendata5))])
#data_normalization
#fifteendata5 <- data.frame(lapply(fifteendata5, function(X) X/X[1]))

#[ order(as.numeric(row.names(data4))), ]
#plot a heatmap
#pheatmap(1/-log10(t(fifteendata5[,2:(ncol(fifteendata5)-1)])), cluster_rows = T, clustering_method = 'ward.D', cluster_cols = F)


fifteendata5_m <- melt(t(fifteendata5[,2:(ncol(fifteendata5))]), id.vars = as.numeric(fifteendata5$time))
fifteensd_m <- melt(t(fifteenstd_df[,2:(ncol(fifteenstd_df))]), id.vars = as.numeric(fifteenstd_df$time))
fifteendata5_m$sd <- fifteensd_m$value

pwarm <- ggplot(data=fifteendata5_m, aes(x=Var2, y=value, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point(aes(color=Var1))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=Var1),  width=.1, position=position_dodge(0.05))+
  facet_grid(Var1 ~ ., scales="free")+
  theme_classic()
pwarm


fifteendata5_m$temperature <- rep("15", nrow(fifteendata5_m))
zerodata5_m$temperature <- rep("0", nrow(zerodata5_m))

data015 <- rbind(fifteendata5_m, zerodata5_m)


pall <- ggplot(data=data015, aes(x=Var2, y=log(value+1), group=temperature)) +
  geom_point(aes(color=temperature))+
  facet_wrap(Var1 ~ ., scales="free")+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=temperature),  width=.1, position=position_dodge(0.05))+
  # scale_color_manual(values=c('green','red'))+
  geom_smooth(method="lm",formula = y ~ poly(x, 3),se = T, aes(color=temperature)) +  facet_wrap(Var1 ~ ., scales="free")+
  scale_color_manual(values=c('steelblue4','orangered3'))+
  theme_classic()
pall



box_plot_conc<-ggplot(data015, aes(x=Var1, y=log10(value+1), fill=temperature)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=14), legend.position="none")+
  geom_boxplot(lwd=.2)+
  scale_fill_manual(values=c('#6699CC','#CC6677'))+
  ylab("log10(concentration+1)")
box_plot_conc




box_plot_conc_PEP<-ggplot(PEP_data, aes(x=Var1, y=log10(value+1), fill=temperature)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=14), legend.position="none")+
  geom_boxplot(lwd=.2)+
  scale_fill_manual(values=c('#6699CC','#CC6677'))+
  ylab("log10(concentration+1)")
box_plot_conc_PEP

ggplot(PEP_data, aes(x=as.character(Var2), y=value, fill=temperature)) + 
  geom_boxplot()


colnames(zerodata5_m) <- c("zeroVar1", "zeroVar2", "zerovalue", "zerotemperature")

# plot the concentration just for timepoint 1
FC_T1 <- fifteendata5[1,]/zerodata5[1,]
FC_T1_m <- melt(sort(FC_T1))

plot(FC_T1_m$value)
plot(density(FC_T1_m$value))

FC_T1_m <- FC_T1_m[-13,]
FC_T1_plot <- ggplot(data=FC_T1_m, aes( x= variable, y=log2(value))) +
  geom_col(size=1, color= "white")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=18))+
  theme(axis.text.y = element_text( hjust=1, size=18))+
  xlab("")+
  ylab("FC")+
  theme(axis.title = element_text(size = 20))  +
  coord_flip()
  
FC_T1_plot

#comparison with model data

#import model data obtained with REMI
predicted_Pflux <- NULL
predicted_Pflux$fluxes <- as.data.frame(read.table("/Users/marcofondi/Dropbox/PhTAC125/AnTReN/Model_integration/REMI_Tania//predicted_PFLux_constrained.txt", sep = '\t', header = T, dec = ','))
predicted_Pflux$mets <- read.table("/Users/marcofondi/Dropbox/PhTAC125/AnTReN/Model_integration/REMI_Tania/mets_for_total_production-flux.txt", sep = '\t', header = F, dec = ',')
predicted_Pflux <- as.data.frame(predicted_Pflux)
met_of_predicted_Pflux <- read.table("/Users/marcofondi/Dropbox/PhTAC125/AnTReN/Model_integration/DeltaFBA-master/lista_me.txt", sep = '\t', header = F, dec = ',')
predicted_Pflux$met <- met_of_predicted_Pflux$V2
predicted_Pflux_sorted <- predicted_Pflux[order(as.numeric(as.character(predicted_Pflux$Var1))),]
FC_T1_m_sorted <- FC_T1_m[order(as.character(FC_T1_m$variable)),]
FC_T1_m_sorted$met <- FC_T1_m_sorted$variable

predicted_Pflux$met <- reorder(predicted_Pflux$met, as.numeric(as.character(predicted_Pflux$Var1)))

#Plot log2FC for time point 1
FC_T1_plot_model <- ggplot(data=predicted_Pflux, aes( x=met, y=log2(as.numeric(as.character(Var1))))) +
  geom_col(size=1, color= "white")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=18))+
  theme(axis.text.y = element_text( hjust=1, size=18))+
  xlab("")+
  ylab("FC")+
  theme(axis.title = element_text(size = 20))  +
  coord_flip()

FC_T1_plot_model

total <- merge(predicted_Pflux,FC_T1_m_sorted,by="met")
total
totalred <- total[c(1,22, 26, 29, 27, 19, 7, 34),]

FC_T1_plot_model_red <- ggplot(data=totalred, aes( x=met, y=log2(as.numeric(as.character(Var1))))) +
  geom_col(size=1, color= "white", fill="steelblue")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=18))+
  theme(axis.text.y = element_text( hjust=1, size=18))+
  xlab("")+
  ylab("FC")+
  theme(axis.title = element_text(size = 20))  +
  coord_flip()

FC_T1_plot_model_red

#compute correlation between model and experimental outliers
corr.test(as.numeric(totalred$value), as.numeric(as.character(totalred$Var1)), method = "spearman")

#Scatterplot of model and experimental outliers
plot(log2(as.numeric(totalred$value)), log(as.numeric(as.character(totalred$Var1))), pch = 21,
     bg = "steelblue",   # Fill color
     col = "white", # Border color
     cex = 2,
     xlab = "Measured concentration (Log2FC)",
     ylab = "Simulated concentration (Log2FC)")     


#compute SMC index for model and experimental data
predicterBin <- as.vector(as.numeric(as.character(total$Var1))>=.95)
measuredBin <- as.vector(as.numeric(as.character(total$value))>=.95)
total$predicterBin <- as.numeric(predicterBin)
total$measuredBin <- as.numeric(measuredBin)
total$smc <- as.numeric(total$predicterBin==total$measuredBin)
sum(total$smc/nrow(total))




# compute FC at T1
FCT1 <- data015_T1$value[1:34]/data015_T1$value[35:nrow(data015_T1)]
MetabolitesFCT1 <- data015_T1$Var1[1:34]

#write concentration out
write.table(data015_T1, file = '/Users/marcofondi/Dropbox/PhTAC125/AnTReN/Model_integration/DeltaFBA-master/Raw_concentrations.txt', sep = '\t',  dec = ',' , row.names = F, col.names = F, quote = F)


#Correlations among metabolite concetrations across time and conditions


data015scatter <- cbind(fifteendata5_m, zerodata5_m)

#plot all possible correlations
pheatmap(cor(zerodata5[,2:(ncol(zerodata5))], fifteendata5[,2:(ncol(fifteendata5))]), cluster_rows = T, cluster_cols = T)
pheatmap(cor(zerodata5[,2:(ncol(zerodata5))], fifteendata5[,2:(ncol(fifteendata5))]))
pheatmap(cor(fifteendata5[,2:(ncol(fifteendata5))],zerodata5[,2:(ncol(zerodata5))] ))
pheatmap(cor(fifteendata5[,2:(ncol(fifteendata5))], fifteendata5[,2:(ncol(fifteendata5))]))
cor_zero_matrix <- cor(zerodata5[,2:(ncol(zerodata5))],zerodata5[,2:(ncol(zerodata5))] )
cor_fifteen_matrix <- cor(fifteendata5[,2:(ncol(fifteendata5))],fifteendata5[,2:(ncol(fifteendata5))] )
corr_differences <- cor_zero_matrix + cor_fifteen_matrix
pheatmap((corr_differences), clustering_method = "ward.D")



pallscatter <- ggplot(data=data015scatter, aes(x=value, y=zerovalue)) +
  geom_point(aes(color=Var1))+
  facet_wrap(Var1 ~ ., scales="free")+
 # scale_color_manual(values=c('green','red'))+
  geom_smooth(method="lm") +
  theme_classic()
pallscatter


## You are now entering the normalized concentration kingdom

# plot normalized concentration 
#fifteendata5 <- data.frame(lapply(fifteendata5, function(X) X/X[1]))

#[ order(as.numeric(row.names(data4))), ]
#plot a heatmap
#pheatmap(1/-log10(t(fifteendata5[,2:(ncol(fifteendata5)-1)])), cluster_rows = T, clustering_method = 'ward.D', cluster_cols = F)


normfifteendata5 <- data.frame(lapply(fifteendata5, function(X) X/X[1]))
normfifteendata5_m <- melt(t(normfifteendata5[,2:(ncol(fifteendata5))]), id.vars = as.numeric(normfifteendata5$time))

#data_normalization
normzerodata5 <- data.frame(lapply(zerodata5, function(X) X/X[1]))
normzerodata5_m <- melt(t(normzerodata5[,2:(ncol(normzerodata5))]), id.vars = as.numeric(normzerodata5$time))


normfifteendata5_m$temperature <- rep("15", nrow(normfifteendata5_m))
normzerodata5_m$temperature <- rep("0", nrow(normzerodata5_m))
data015_norm <- rbind(normfifteendata5_m, normzerodata5_m)


# Figure 1A
pall_norm <- ggplot(data=data015_norm, aes(x=Var2, y=value, group=temperature)) +
  geom_point(aes(color=temperature))+
  geom_smooth(method="lm",formula = y ~ poly(x, 4),se = T, aes(color=temperature)) +  facet_wrap(Var1 ~ ., scales="free")+
  scale_color_manual(values=c('#6699CC','#CC6677'))+
  theme_classic()

pall_norm <- pall_norm + labs(x=NULL, y=NULL)
pall_norm <- pall_norm + theme_bw(base_size=10)
pall_norm <- pall_norm + theme(strip.background=element_blank())
pall_norm <- pall_norm + theme(axis.text.x=element_text(angle=90, vjust=0.5))
pall_norm <- pall_norm + theme(panel.grid.major.x=element_blank())
pall_norm <- pall_norm + theme(panel.grid.major.y=element_blank())
pall_norm <- pall_norm + theme(panel.grid.minor.y=element_blank())
pall_norm



#correlation in normalized data
colnames(normzerodata5_m) <- c("zeroVar1", "zeroVar2", "zerovalue", "zerotemperature")
data015scatter_norm <- cbind(normfifteendata5_m, normzerodata5_m)
pallscatter <- ggplot(data=data015scatter_norm, aes(x=value, y=zerovalue)) +
  geom_point(aes(color=Var1))+
  facet_wrap(Var1 ~ ., scales="free")+
  # scale_color_manual(values=c('green','red'))+
  geom_smooth(method="lm") +
  theme_classic()
pallscatter

#plot all possible  correlations of normalized data
cor_normzero_matrix <- cor(normzerodata5[,2:(ncol(normzerodata5))],normzerodata5[,2:(ncol(normzerodata5))] )
cor_normfifteen_matrix <- cor(normfifteendata5[,2:(ncol(normfifteendata5))],normfifteendata5[,2:(ncol(normfifteendata5))] )
corr_norm_differences <- cor_normzero_matrix + cor_normfifteen_matrix
pheatmap((corr_norm_differences), clustering_method = "ward.D")

correlation_0vs15 <- cor(zerodata5[,2:(ncol(zerodata5))], fifteendata5[,2:(ncol(fifteendata5))])
diagonal_correlation <- diag(correlation_0vs15)
dfm <- melt(diagonal_correlation)

pheatmap((correlation_0vs15))

plot(sort(diagonal_correlation), xaxt = "n", xlab='' , ylab = 'Pearson correlation coeffcient',pch  = 1,
     col = "blue")          # Box color

axis(1, at=1:length(sort(diagonal_correlation)),  labels=names(sort(diagonal_correlation)), las=2)


#ratio between temperature 0 divided by 15

ratio_0_15 <- zerodata5/fifteendata5

ratio_0_15_m <- melt(t(ratio_0_15[,2:(ncol(ratio_0_15))]), id.vars = as.numeric(ratio_0_15$OD))



pall_ratio <- ggplot(data=ratio_0_15_m, aes(x=Var2, y=value, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point(aes(color=Var1))+
  facet_wrap(Var1 ~ ., scale='free')+
  #scale_y_continuous(limits = c(0, 3))+
  theme_classic()
pall_ratio

#correlation with normalized data

correlation_0vs15_norm <- cor(normzerodata5[,2:(ncol(normzerodata5))], normfifteendata5[,2:(ncol(normfifteendata5))], method = "pearson")
diagonal_correlation_norm <- diag(correlation_0vs15_norm)

col<- colorRampPalette(c("red", "white", "blue"))(20)
corrplot(correlation_0vs15_norm, order="hclust", col = COL2('RdYlBu'), method = 'square', type="lower")

plot(sort(diagonal_correlation_norm), xaxt = "n", xlab='' , cex=2, ylab = 'Pearson correlation coeffcient',pch  = 16, size=3,col = "black")          # Box color
axis(1, at=1:length(sort(diagonal_correlation)),  labels=names(sort(diagonal_correlation)), las=2)
diagonal_correlation_norm_m <- melt(diagonal_correlation_norm)
diagonal_correlation_norm_m$metabolite <- rownames(diagonal_correlation_norm_m)
diagonal_correlation_norm_m <- diagonal_correlation_norm_m[order(diagonal_correlation_norm_m$value),]
diagonal_correlation_norm_m$metabolite <- factor(diagonal_correlation_norm_m$metabolite, levels = diagonal_correlation_norm_m$metabolite[order(diagonal_correlation_norm_m$value)])

correlation_plot <- ggplot(data=diagonal_correlation_norm_m, aes( x= metabolite, y=value)) +
  geom_point(size=4, color= "#999933")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=14))+
  xlab("")+
  ylab("Correlation (PPM)")

correlation_plot



cor.test(normzerodata5$PEP, normfifteendata5$PEP)
library(heatmaply)


cor.test.p <- function(x,w){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], w[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(w))
  z
}


  p <- cor.test.p(normzerodata5[,2:(ncol(normzerodata5)-1)], normfifteendata5[,2:(ncol(normfifteendata5)-1)])
  
  heatmaply_cor(
    correlation_0vs15_norm,
    node_type = "scatter",
    point_size_mat = -log(p), 
    point_size_name = "-log10(p-value)",
    label_names = c("x", "y", "Correlation")
  )

#correlation within temp
correlation_0vs0_norm <- cor(normzerodata5[,2:(ncol(normzerodata5))], normzerodata5[,2:(ncol(normzerodata5))], method = "pearson")
pheatmap(correlation_0vs0_norm)
corrplot(correlation_0vs0_norm, type="upper", order="hclust", col=col)

correlation_15vs15_norm <- cor(normfifteendata5[,2:(ncol(normfifteendata5))], normfifteendata5[,2:(ncol(normfifteendata5))], method = "pearson")
pheatmap(correlation_15vs15_norm)
corrplot(correlation_15vs15_norm, type="lower", order="hclust", col=col)

#test whether correlation matrices (ovs15, 0vs0 and 15vs15) are statistically different
cor.test(correlation_15vs15_norm, correlation_0vs15_norm)
t.test(correlation_15vs15_norm, correlation_0vs15_norm, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)

grid.newpage()
#grid.arrange(GG_uptake_0_plot,GG_uptake_15_plot, ncol=2, nrow=1)
grid.arrange( correlation_plot,box_plot_conc, pall_norm, ncol=2, nrow=2 )


