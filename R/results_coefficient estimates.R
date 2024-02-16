
library(jagsUI)
library(ggplot2)

out <- readRDS(here::here("outputs","resu_modelCovars3Pool.rds"))

out$summary

# Scale of effect
table(out$sims.list$scale_for) / out$mcmc.info$n.samples

# Parameters coefficient estimates -----------------------------------------------------
coefs.tab <- as.data.frame(out$summary[1:(grep("w",rownames(out$summary))[1] - 1),c("mean","sd","2.5%","97.5%")])

coefs.tab$par <- as.vector(do.call(rbind,strsplit(rownames(coefs.tab),"[",fixed=T))[,1])

coefs.tab$sp <- NA
coefs.tab[grep("]",rownames(coefs.tab)),"sp"] <- dimnames(Y.DA)[[3]]
coefs.tab[-grep("]",rownames(coefs.tab)),"sp"] <- "M-L mammals"


# psi ~ forest1000 --------------------------------------------------------
# get parameter estimates
psi.for1000 <- data.frame(
  sp=c("M-L Mammals",mammInfo$sp.id),
  mean=c(out$mean$mu.delta1,out$mean$delta1),
  lcl=c(out$q2.5$mu.delta1,out$q2.5$delta1),
  ucl=c(out$q97.5$mu.delta1,out$q97.5$delta1),
  signif=!c(out$overlap0$mu.delta1,out$overlap0$delta1)
)

# keep only detected species
psi.for1000 <- psi.for1000[which(psi.for1000$sp %in% names(which(colSums(ncap)>0)) | psi.for1000$sp=="M-L Mammals"),]

psi.for1000 <- psi.for1000[order(psi.for1000$mean,decreasing=T),]
psi.for1000$sp <- factor(psi.for1000$sp, levels=c(rev(as.character(psi.for1000$sp))[-9],"M-L Mammals"))

fig.psiFor1000 <- ggplot(data=psi.for1000,aes(y=mean,x=sp,ymin=lcl,ymax=ucl,col=signif)) +
  geom_point(size=ifelse(psi.for1000$sp=="M-L Mammals",5,3),shape=ifelse(psi.for1000$sp=="M-L Mammals",17,19)) + 
  geom_errorbar(size=ifelse(psi.for1000$sp=="M-L Mammals",1.4,1),width=0.5) +
  geom_hline(yintercept=0,col="gray60",linetype="dashed") +
  #  geom_hline(yintercept=c(out$mean$mu.delta1,out$q2.5$mu.delta1,out$q97.5$mu.delta1),
  #             col="red",linetype=c("solid","dashed","dashed"),size=1) +
  scale_color_manual(labels=c("non-significant","significant"),values=c("gray60","blue")) +
  labs(x="",y="Effect of % NF forests (1km buffer)\n in the carcass suitability",color="") +
  theme_classic() + coord_flip() +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),
        legend.position="top",legend.text=element_text(size=14))

ggsave(file.path("figs","Fig.Psi_Coefs_Forest.png"),fig.psiFor1000,
       width=16,height=20,units="cm",dpi=320)
