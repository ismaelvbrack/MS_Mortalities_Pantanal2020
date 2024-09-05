
library(jagsUI)
library(ggplot2)
library(ggrepel)
source("https://raw.githubusercontent.com/ismaelvbrack/predictJAGS_func/main/func_predictJAGS.R")


out <- readRDS(here::here("outputs","resu_modelCovars3Pool.rds"))

# Carcass occurrence ~ forest (1000m buffer) ------------------------------

new.for1000 <- data.frame(forest1000=seq(min(plotsCovs.s[,"forest1000"]),max(plotsCovs.s[,"forest1000"]),,100))

pred.thetaFor.multi <- predictJAGS(fm=out,params=c("mu.delta0","mu.delta1"),
                                   newdata=new.for1000,link="logit",quants=c(0.025,0.5,0.975),appendData=T)

pred.thetaFor.multi$forest1000 <- new.for1000[,"forest1000"]*
  attr(plotsCovs.s, "scaled:scale")["forest1000"]+
  attr(plotsCovs.s, "scaled:center")["forest1000"]

names(pred.thetaFor.multi)[2:4] <- c("lower","median","upper")


pred.thetaFor.sp <- matrix(NA,ncol=nspp,nrow=nrow(new.for1000),
                           dimnames=list(NULL,mammInfo$sp.id))
for(i in 1:nspp){
  pred.thetaFor.sp[,i] <- plogis(out$mean$delta0[i] + (out$mean$delta1[i]*new.for1000$forest1000))
}

pred.thetaFor.sp <- data.frame(
  taxa=rep(mammInfo$sp.id,  
           each=nrow(new.for1000)),
  forest1000=rep(new.for1000$forest1000,nspp),
  est=as.vector(pred.thetaFor.sp)
)
pred.thetaFor.sp$forest1000 <- pred.thetaFor.sp$forest1000*
  attr(plotsCovs.s, "scaled:scale")["forest1000"]+
  attr(plotsCovs.s, "scaled:center")["forest1000"]

# exclude not-detected species
pred.thetaFor.sp <- pred.thetaFor.sp[pred.thetaFor.sp$taxa %in% names(which(colSums(ncap)>0)),]


#* Figure!
signifs <- data.frame(est1=tapply(pred.thetaFor.sp$est,pred.thetaFor.sp$taxa,function(x) x[100]),
           taxa=unique(pred.thetaFor.sp$taxa),
           signif=!out$overlap0$delta1[colSums(ncap)>0],
           cols="gray")
signifs[which(signifs$taxa=="Bladic"),"signif"] <- FALSE

signifs[which(signifs$signif==1),"cols"] <- rep(c("red3","orange"),ceiling(sum(signifs$signif==1)/2))
signifs[which(signifs$taxa=="Nasnas"),"cols"] <- "red3"

fig.psi.for1000 <- 
  ggplot(data=pred.thetaFor.sp,aes(x=forest1000,y=est,group=taxa,label=taxa,col=taxa)) + 
  geom_line(size=1,alpha=0.6) +
  theme(axis.title=element_text(size=18),axis.text=element_text(size=16),plot.margin=margin(20,80,20,20)) + # plot.margin=margin(20,40,20,20)
  theme_classic() + #coord_cartesian(ylim=c(0,100)) +
  labs(x="% of non-flooded forests (1km buffer)",y="Prob. of carcass suitability (\U03A8)") +
  scale_x_continuous(breaks=seq(.1,.9,0.2),labels=seq(10,90,20),expand=c(0,0.1)) +
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
  scale_color_manual(values=signifs$cols) +
  geom_line(data=pred.thetaFor.multi, aes(x=forest1000,y=median),
            size=1.6,color="darkgreen",inherit.aes=F) +
  geom_ribbon(data=pred.thetaFor.multi, aes(x=forest1000,y=median,ymax=upper,ymin=lower),
              alpha=0.3, fill="darkgreen",inherit.aes=F) +
  geom_text_repel(data=signifs,
                  aes(x=rep(1,sum(colSums(ncap)>0)),y=est1,label=taxa),inherit.aes=F,
                  fontface=ifelse(signifs$signif==FALSE,"plain","bold"),
                  nudge_x=.2,direction="y",hjust="right",
                  col=signifs$cols,#ifelse(out$overlap0$delta1[colSums(ncap)>0],"gray50","black"),
                  segment.size=.3,segment.alpha=.7,size=5,
                  box.padding=.4) +
  theme(legend.position="none",axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))


 # ggplot(data=pred.thetaFor.sp,aes(x=forest1000,y=est,group=taxa,label=taxa)) + 
 #  geom_line(col="gray40",size=1,alpha=0.4) +
 #  theme(axis.title=element_text(size=18),axis.text=element_text(size=16),plot.margin=margin(20,80,20,20)) + # plot.margin=margin(20,40,20,20)
 #  theme_classic() + #coord_cartesian(ylim=c(0,100)) +
 #  labs(x="% of non-flooded forests (1km buffer)",y="Prob. of carcass suitability (\U03A8)") +
 #  scale_x_continuous(breaks=seq(.1,.9,0.2),labels=seq(10,90,20),expand=c(0,0.1)) +
 #  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2)) +
 #  geom_line(data=pred.thetaFor.multi, aes(x=forest1000,y=median),
 #            size=1.6,color="darkgreen",inherit.aes=F) +
 #  geom_ribbon(data=pred.thetaFor.multi, aes(x=forest1000,y=median,ymax=upper,ymin=lower),
 #              alpha=0.4, fill="darkgreen",inherit.aes=F) +
 #  geom_text_repel(data=,
 #  aes(x=rep(1,sum(colSums(ncap)>0)),y=est1,label=taxa),inherit.aes=F,
 #  fontface=ifelse(out$overlap0$delta1[colSums(ncap)>0],"plain","bold"),
 #  nudge_x=.2,direction="y",hjust="right",
 #  col=ifelse(out$overlap0$delta1[colSums(ncap)>0],"gray50","black"),
 #  segment.size=.3,segment.alpha=.7,size=4,
 #  box.padding=.4) +
 #  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"))

ggsave(here::here("figs","Fig.Psi-Forest.png"),fig.psi.for1000,
       width=18,height=16,units="cm",dpi=320)



# Carcass abundance ~ deltaBurn + tanque ----------------------------------

new.dBurnT <- data.frame(deltaBurn=rep(seq(min(plotsCovs.s[,"deltaBurn"]),max(plotsCovs.s[,"deltaBurn"]),,100),2),
                         tanque=rep(0:1,each=100))

#out$sims.list$mean.alpha0 <- apply(out$sims.list$mu.alpha0,1,mean)
pred.lambda.multi <- predictJAGS(fm=out,params=c("mu.alpha0","alpha1","alpha2"),
                                 newdata=new.dBurnT,link="log",quants=c(0.025,0.5,0.975),appendData=T)

pred.lambda.multi$deltaBurn <- new.dBurnT[,"deltaBurn"]*
  attr(plotsCovs.s, "scaled:scale")["deltaBurn"]+
  attr(plotsCovs.s, "scaled:center")["deltaBurn"]

pred.lambda.multi$tanque <- factor(pred.lambda.multi$tanque,labels=c("no-tanque","tanque"))

names(pred.lambda.multi)[2:4] <- c("lower","median","upper")


pred.lambda.sp <- matrix(NA,ncol=nspp,nrow=nrow(new.dBurnT),
                         dimnames=list(NULL,mammInfo$sp.id))
for(i in 1:nspp){
  pred.lambda.sp[,i] <- exp(out$mean$alpha0[i] + (out$mean$alpha1*new.dBurnT$deltaBurn) + (out$mean$alpha2*new.dBurnT$tanque))
}

pred.lambda.sp <- data.frame(
  taxa=rep(mammInfo$sp.id,  
           each=nrow(new.dBurnT)),
  deltaBurn=rep(new.dBurnT$deltaBurn,nspp),
  tanque=rep(c("no-tanque","tanque"),each=100),
  est=as.vector(pred.lambda.sp)
)
pred.lambda.sp$deltaBurn <- pred.lambda.sp$deltaBurn*
  attr(plotsCovs.s, "scaled:scale")["deltaBurn"]+
  attr(plotsCovs.s, "scaled:center")["deltaBurn"]

#* Figure!
fig.dburn.tanque <- ggplot(data=pred.lambda.sp,aes(x=deltaBurn,y=est,group=interaction(taxa,tanque),col=as.factor(tanque),label=taxa)) + 
  #geom_line(size=1,alpha=0.2) +
  theme(axis.title=element_text(size=18),axis.text=element_text(size=16),plot.margin=margin(20,80,20,20)) + # plot.margin=margin(20,40,20,20)
  theme_classic() + #coord_cartesian(ylim=c(0,100)) +
  labs(x="Wildfire severity (ΔNBR)",y="Carcasses / ha",color="Tanque") +
  scale_color_manual(values=c("brown1","cyan3"),guide=guide_legend(override.aes=list(alpha=1,size=1.4,color=c("red4","blue4")))) +
  scale_x_continuous(expand=c(0,0.1)) +
  #scale_y_continuous(breaks=seq(0,2.5,0.5)) +
  #coord_cartesian(ylim=c(0,2.2)) + 
  facet_wrap(vars(tanque),nrow=2,scales="free",labeller=as_labeller(c(`no-tanque`="without AWB",`tanque`="with AWB"))) +
  geom_line(data=pred.lambda.multi[which(pred.lambda.multi$tanque=="no-tanque"),], aes(x=deltaBurn,y=mean),
            size=1.6,inherit.aes=F,color="red4") + 
  geom_line(data=pred.lambda.multi[which(pred.lambda.multi$tanque=="tanque"),], aes(x=deltaBurn,y=mean),
            size=1.6,inherit.aes=F,color="blue4") +
  geom_ribbon(data=pred.lambda.multi[which(pred.lambda.multi$tanque=="no-tanque"),], aes(x=deltaBurn,y=mean,ymax=upper,ymin=lower),
              alpha=0.3, fill="red4",inherit.aes=F) +
  geom_ribbon(data=pred.lambda.multi[which(pred.lambda.multi$tanque=="tanque"),], aes(x=deltaBurn,ymax=upper,ymin=lower),
              alpha=0.3, fill="blue4",inherit.aes=F) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16,face="bold"),
        strip.text=element_text(size=14,face="bold"))

ggsave(here::here("figs","Fig.Lambda-dBurnAWB.png"),fig.dburn.tanque,
       width=16,height=24,units="cm",dpi=320)
