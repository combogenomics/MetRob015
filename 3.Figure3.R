#caricamento librerie
library(ggplot2)
require(pheatmap)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ggplot2)
library(reshape2)
library(heatmap3)
library(forcats)
library(grid)
library(gridExtra)
library(corrplot)
library(readxl)
setwd("/Users/marcofondi//Dropbox//PhTAC125/AnTReN/metabolomics/")
#importare file
met_ext_0 <- read.table("met_ext_0 _not0.txt", sep = '\t', header = T, dec = ',')
met_ext_15 <- read.table("met_ext_15_not0.txt", sep = '\t', header = T, dec = ',')

#Rimozione OD
met_ext_0 = met_ext_0[, -c(1)]
met_ext_15 = met_ext_15[, -c(1)]
#calcolo medie
#media dei metaboliti a 15?C
met_ext_15_media <- do.call(rbind,
                            lapply(seq(1, nrow(met_ext_15), 2), function(i){
                              x <- met_ext_15[ i:(i + 1), , drop = FALSE]
                              res <- rbind( colSums(x)/2)
                              rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                              res
                            }))

#media dei metaboliti a 0?C
met_ext_0_media <- do.call(rbind,
                           lapply(seq(1, nrow(met_ext_0), 2), function(i){
                             x <- met_ext_0[ i:(i + 1), , drop = FALSE]
                             res <- rbind( colSums(x)/2)
                             rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                             res
                           }))
#conversione formato in dataset
met_ext_15_media2 <- data.frame(met_ext_15_media)
rownames(met_ext_15_media2)<- c(1:5)

met_ext_0_media2 <- data.frame(met_ext_0_media)
rownames(met_ext_0_media2)<-c(1:5)
#calcolo deviazione standard
met_ext_15_std <- do.call(rbind,
                          lapply(seq(1, nrow(met_ext_15), 2), function(i){
                            x <- met_ext_15[ i:(i + 1), , drop = FALSE]
                            res <- rbind( apply(x, 2, sd, na.rm = TRUE))
                            rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                            res
                          }))

met_ext_0_std <- do.call(rbind,
                         lapply(seq(1, nrow(met_ext_0), 2), function(i){
                           x <- met_ext_0[ i:(i + 1), , drop = FALSE]
                           res <- rbind( apply(x, 2, sd, na.rm = TRUE))
                           rownames(res)[ nrow(res) ] <- paste(rownames(x), collapse = "_")
                           res
                         }))
#conversione formato in dataset
met_ext_15_std2 <- data.frame(met_ext_15_std)
row.names(met_ext_15_std2) <- c(1:5)

met_ext_0_std2 <- data.frame(met_ext_0_std)
row.names(met_ext_0_std2) <- c(1:5)


met_ext_15_media2_m <- melt(t(met_ext_15_media2[,1:(ncol(met_ext_15_media2))]), id.vars = as.numeric(met_ext_15_media2$time))
met_ext_15_std2_m <- melt(t(met_ext_15_std2[,1:(ncol(met_ext_15_std2))]), id.vars = as.numeric(met_ext_15_std2$time))
met_ext_15_media2_m$sd <- met_ext_15_std2_m$value


met_ext_0_media2_m <- melt(t(met_ext_0_media2[,1:(ncol(met_ext_0_std2))]), id.vars = as.numeric(met_ext_0_std2$time))
met_ext_0_std2_m <- melt(t(met_ext_0_std2[,1:(ncol(met_ext_0_std2))]), id.vars = as.numeric(met_ext_0_std2$time))
met_ext_0_media2_m$sd <- met_ext_0_std2_m$value


met_ext_15_media2_m$temperature <- rep("15", nrow(met_ext_15_media2_m))
met_ext_0_media2_m$temperature <- rep("0", nrow(met_ext_0_media2_m))

data015_extr <- rbind(met_ext_15_media2_m, met_ext_0_media2_m)

pall_extra <- ggplot(data=data015_extr, aes(x=Var2, y=value, group=temperature)) +
  geom_line(aes(color=temperature))+
  geom_point(aes(color=temperature))+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=temperature),  width=.1, position=position_dodge(0.05))+
  facet_wrap(Var1 ~ ., scales="free")+
  scale_color_manual(values=c('green','red'))+
  theme_classic()
pall_extra

#rownames(data015_extr) <- data015_extr$Var1

#get all the metabolites that are common to 0 and 15

shared_mets <- intersect(colnames(met_ext_0), colnames(met_ext_15))
# met_ext_0[, colnames(met_ext_0) == shared_mets]
# data015_extr[data015_extr$Var1 ==  as.factor(shared_mets),]

data015_extr_shared <- data015_extr[data015_extr$Var1 %in% shared_mets, ]  


pall_extra_shared <- ggplot(data=data015_extr_shared, aes(x=Var2, y=value, group=temperature)) +
  geom_point(aes(color=temperature))+
  facet_wrap(Var1 ~ ., scales="free")+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=temperature),  width=.1, position=position_dodge(0.05))+
  # scale_color_manual(values=c('green','red'))+
  geom_smooth(method="lm",formula = y ~ poly(x, 2),se = TRUE, aes(color=temperature)) +  facet_wrap(Var1 ~ ., scales="free")+
  scale_color_manual(values=c('#6699CC','#CC6677'))+
  theme_classic()
pall_extra_shared

pall_extra_shared <- pall_extra_shared + labs(x=NULL, y=NULL)
pall_extra_shared <- pall_extra_shared + theme_bw(base_size=10)
pall_extra_shared <- pall_extra_shared + theme(strip.background=element_blank())
pall_extra_shared <- pall_extra_shared + theme(axis.text.x=element_text(angle=90, vjust=0.5))
pall_extra_shared <- pall_extra_shared + theme(panel.grid.major.x=element_blank())
pall_extra_shared <- pall_extra_shared + theme(panel.grid.major.y=element_blank())
pall_extra_shared <- pall_extra_shared + theme(panel.grid.minor.y=element_blank())
pall_extra_shared

box_plot_conc_extra<-ggplot(data015_extr_shared, aes(x=Var1, y=log10(value), fill=temperature)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=14), legend.position="none")+
  geom_boxplot(lwd=.2)+
  scale_fill_manual(values=c('#6699CC','#CC6677'))+
  ylab("log10(concentration)")
box_plot_conc_extra

data015scatter <- cbind(met_ext_0[, shared_mets], met_ext_15[, shared_mets])

met_ext_0_shared <- met_ext_0_media2[, shared_mets]
met_ext_15_shared <- met_ext_15_media2[, shared_mets]

correlation_0vs15_extra <- cor(met_ext_0_shared[,1:ncol(met_ext_0_shared)], met_ext_15_shared[,1:ncol(met_ext_15_shared)])
diagonal_correlation_extra <- diag(correlation_0vs15_extra)

pheatmap((correlation_0vs15_extra))

plot(sort(diagonal_correlation_extra), xaxt = "n", xlab='' , ylab = 'Pearson correlation coeffcient',pch  = 16,
     col = "black", cex=2)          # Box color
axis(1, at=1:length(sort(diagonal_correlation_extra)),  labels=names(sort(diagonal_correlation_extra)), las=2)
cor.test(met_ext_0_shared$X2oxoglutarate, met_ext_15_shared$X2oxoglutarate)

diagonal_correlation_extra_m <- melt(diagonal_correlation_extra)
diagonal_correlation_extra_m$metabolite <- rownames(diagonal_correlation_extra_m)
diagonal_correlation_extra_m <- diagonal_correlation_extra_m[order(diagonal_correlation_extra_m$value),]
diagonal_correlation_extra_m$metabolite <- factor(diagonal_correlation_extra_m$metabolite, levels = diagonal_correlation_extra_m$metabolite[order(diagonal_correlation_extra_m$value)])

correlation_plot_extra <- ggplot(data=diagonal_correlation_extra_m, aes( x= metabolite, y=value)) +
  geom_point(size=4, color= "#999933")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=14))+
  xlab("")+
  ylab("Correlation (PPM)")
correlation_plot_extra
#normalized data by first measure

met_ext_0_norm <- data.frame(lapply((met_ext_0_shared+0.1), function(X) X/(X[1])))
#met_ext_0_norm$time <- rownames(met_ext_0_norm)
met_ext_0_norm_m <- melt(t(met_ext_0_norm), id.vars = as.numeric(rownames(met_ext_0_norm$time)))
met_ext_0_norm_m$temperature<- rep("0", nrow(met_ext_0_norm_m))

met_ext_15_norm <- data.frame(lapply((met_ext_15_shared+0.1), function(X) X/(X[1])))
#$met_ext_15_norm$time <- rownames(met_ext_15_norm)
met_ext_15_norm_m <- melt(t(met_ext_15_norm), id.vars = as.numeric(rownames(met_ext_15_norm)))
met_ext_15_norm_m$temperature<- rep("15", nrow(met_ext_15_norm_m))

all_extra_norma <- rbind(met_ext_0_norm_m, met_ext_15_norm_m)

pall_norm_extra <- ggplot(data=all_extra_norma, aes(x=Var2, y=value, group=temperature)) +
  geom_point(aes(color=temperature))+
    #geom_errorbar(aes(ymin=value-sd, ymax=value+sd, color=temperature),  width=.1, position=position_dodge(0.05))+
  # scale_color_manual(values=c('green','red'))+
  geom_smooth(method="lm",formula = y ~ poly(x, 4),se = F, aes(color=temperature)) +  facet_wrap(Var1 ~ ., scales="free")+
  scale_color_manual(values=c('#6699CC','#CC6677'))+
  theme_classic()
pall_norm_extra

pall_norm_extra <- pall_norm_extra + labs(x=NULL, y=NULL)
pall_norm_extra <- pall_norm_extra + theme_bw(base_size=10)
pall_norm_extra <- pall_norm_extra + theme(strip.background=element_blank())
pall_norm_extra <- pall_norm_extra + theme(axis.text.x=element_text(angle=90, vjust=0.5))
pall_norm_extra <- pall_norm_extra + theme(panel.grid.major.x=element_blank())
pall_norm_extra <- pall_norm_extra + theme(panel.grid.major.y=element_blank())
pall_norm_extra <- pall_norm_extra + theme(panel.grid.minor.y=element_blank())
pall_norm_extra

ks.test(met_ext_0$Hypoxanthine, met_ext_15$Hypoxanthine, alternative= "greater")


data015_T1_extr <- data015_extr[data015_extr$Var2==1,]

FCT1_ext <- data015_T1_extr$value[1:20]/data015_T1_extr$value[21:nrow(data015_T1_extr)]
MetabolitesFCT1_extr <- data015_T1_extr$Var1[1:20]

FCT1_df_extr <- NULL
ratio_production_flux_extr <- NULL
FCT1_df_extr$meta <- MetabolitesFCT1_extr
FCT1_df_extr$FC <- FCT1_ext
FCT1_df_extr <- as.data.frame(FCT1_df_extr)



