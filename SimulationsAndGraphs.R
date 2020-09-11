setwd("E:/Letterv3/")
library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)
library(svglite)

ChevalCount <- read_excel("RawData.xlsx",
                          sheet = 1, col_names = T)
DiaoCount <- read_excel("RawData.xlsx",
                        sheet = 2, col_names = T)
DiaoLayer <- read_excel("RawData.xlsx",
                        sheet = 3, col_names = T)

gen_wilcox<-function(sd2=1, n=100){
  X <- rnorm(n, 0, 1)
  Y <- rnorm(n, 0, sd2)
  return(wilcox.test(X,Y)$p.value<0.05)
}
gen_boot_H0<-function(sd2=1, n=100, N=1000){
  X <- rnorm(n, 0, 1)
  Y <- rnorm(n, 0, sd2)
  mediandiff<-median(X)-median(Y)
  boots<-replicate(N, median(sample(X,n, replace=T))-median(sample(Y,n, replace=T))-mediandiff)
  above <- sum(abs(boots)>=abs(mediandiff))
  return((above+1)/(N+1)<0.05)
}
gen_wilcox_H0_beta<-function(sd2=1, n=100, N=1000){
  X <- rbeta(n, 1, 3)
  Y <- rnorm(n, 1 - 1/2^(1/3), sqrt(.0375))
  return(wilcox.test(X,Y)$p.value<0.05)
}
gen_boot_H0_beta<-function(sd2=1, n=100, N=1000){
  X <- rbeta(n, 1, 3)
  Y <- rnorm(n, 1 - 1/2^(1/3), sqrt(.0375))
  mediandiff<-median(X)-median(Y)
  boots<-replicate(N, median(sample(X,n, replace=T))-median(sample(Y,n, replace=T))-mediandiff)
  above <- sum(abs(boots)>=abs(mediandiff))
  return((above+1)/(N+1)<0.05)
}
mcp_ci <- function(success, trials, alpha=0.05){
  ## Copyright (C) 2001 Frank E Harrell Jr
  ## Modified by Matthew G Johnston 2020 from binconf in the Hmisc package
  ## Distributed under GNU General Public License v2 (or later) without any warranty. See the GNU General Public License for more details.
  zcrit <-  - qnorm(alpha/2)
  z2 <- zcrit * zcrit
  mc_p <- success/trials
  cl <- (mc_p + z2/2/trials + c(-1, 1) * zcrit *
           sqrt((mc_p * (1 - mc_p) + z2/4/trials)/trials))/(1 + z2/trials)
  if(success == 1)
    cl[1] <-  - log(1 - alpha)/trials
  if(success == (trials - 1))
    cl[2] <- 1 + log(1 - alpha)/trials
  return(cbind(mc_p, lower_ci=cl[1], upper_ci=cl[2]))
}
boot_median<- function(data1, data2, N=10000){
  mediandiff<-median(data1)-median(data2)
  boots<-replicate(N, median(sample(data1,length(data1), replace=T))-median(sample(data2,length(data2),  replace=T))-mediandiff)
  above <- sum(abs(boots)>=abs(mediandiff))
  mcp<-format(round(mcp_ci(above+1,N+1),3),nsmall=3)
  labeltext1<- substitute(atop(paste(hat(italic("p")))~phantom()==phantom()~AAA,"95% CI"~BBB~CCC), list(AAA=mcp[1], BBB =paste0("[",mcp[2],", "), CCC= paste0(mcp[3],"]")))
  plot<-ggplot(data.frame(median=boots), aes(x=abs(median)))+
    geom_histogram(alpha=.7, binwidth = 1, aes(y=..density..))+
    geom_vline(xintercept = abs(mediandiff), linetype=2, colour = "red4")+
    theme_bw()+
    xlab(expression(atop("|"~hat(paste(theta, "*"))~-~hat(theta)~"|",theta==Delta[median])))+
    ylab("Density")+
    annotate("text",x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,label=labeltext1)+
    theme(text=element_text(size=25))
  return(plot)
}
boot_mean<- function(data1, data2, N=10000){
  meandiff<-mean(data1)-mean(data2)
  boots<-replicate(N, mean(sample(data1,length(data1),  replace=T))-mean(sample(data2,length(data2),  replace=T))-meandiff)
  above <- sum(abs(boots)>=abs(meandiff))
  mcp<-format(round(mcp_ci(above+1,N+1),3),nsmall=3)
  labeltext1<- substitute(atop(paste(hat(italic("p")))~phantom()==phantom()~AAA,"95% CI"~BBB~CCC), list(AAA=mcp[1], BBB =paste0("[",mcp[2],", "), CCC= paste0(mcp[3],"]")))
  plot<-ggplot(data.frame(mean=boots), aes(x=abs(mean)))+
    geom_histogram(alpha=.7, bins = 10, aes(y=..density..))+
    geom_vline(xintercept = abs(meandiff), linetype=1, colour = "red4")+
    theme_bw()+
    xlab(expression(atop("|"~hat(paste(theta, "*"))~-~hat(theta)~"|",theta==Delta[mu])))+
    ylab("Density")+
    annotate("text",x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,label=labeltext1)+
    theme(text=element_text(size=25),
          axis.title.x =  element_text(size=25))
  return(plot)
}

N<-1000
seed <- 9

set.seed(seed)
equal<-sum(replicate(N, gen_wilcox(n=100)))
unequal<-sum(replicate(N, gen_wilcox(sd2=5,n=100)))
beta<-sum(replicate(N, gen_wilcox_H0_beta(sd2=5,n=100)))
mcp_ci(equal,N)*100
mcp_ci(unequal,N)*100
mcp_ci(beta,N)*100

set.seed(seed)
equal<-sum(replicate(N, gen_boot_H0(n=100)))
unequal<-sum(replicate(N, gen_boot_H0(sd2=5,n=100)))
beta<-sum(replicate(N, gen_boot_H0_beta(sd2=5,n=100)))
mcp_ci(equal,N)*100
mcp_ci(unequal,N)*100
mcp_ci(beta,N)*100

###Cheval
ChevalCount$treatment<-factor(ChevalCount$treatment, levels=c("Mock","Chitin"))

vline_df <- data.frame(treatment = levels(ChevalCount$treatment),
                       Medians = tapply(X = ChevalCount$movement, INDEX = ChevalCount$treatment,
                                        FUN = median))
vline2_df <- data.frame(treatment = levels(ChevalCount$treatment),
                        Means = tapply(X = ChevalCount$movement, INDEX = ChevalCount$treatment,
                                       FUN = mean))

Figure1<-ggplot(ChevalCount,
                aes(x=movement, y=..density..))+
  geom_histogram(alpha=.7, binwidth = 1)+
  facet_wrap(~treatment, ncol=1)+
  coord_cartesian(xlim=c(0,40))+
  xlab("Number of cells GFP moves")+
  ylab("Density")+theme_bw()+
  geom_rug(aes(y=0),position = position_jitter(height = 0), alpha=0.2)+
  geom_vline(data = vline_df, aes(xintercept = Medians), linetype = 2, colour = "red4")+
  geom_vline(data = vline2_df, aes(xintercept = Means), linetype = 1, colour = "red4")+
  theme(text=element_text(size=25))
data1 <- ChevalCount[ChevalCount$treatment=="Mock",]$movement
data2 <- ChevalCount[ChevalCount$treatment=="Chitin",]$movement
Figure2 <- boot_median(data1,data2)
Figure3 <- boot_mean(data1,data2)

###Diao

DiaoCount$genotype<-factor(DiaoCount$genotype, levels(as.factor(DiaoCount$genotype))[c(2,1)])
vline_df <- data.frame(genotype = levels(DiaoCount$genotype),
                       Medians = tapply(X = DiaoCount$movement, INDEX = DiaoCount$genotype,
                                        FUN = median))
vline2_df <- data.frame(genotype = levels(DiaoCount$genotype),
                        Means = tapply(X = DiaoCount$movement, INDEX = DiaoCount$genotype,
                                       FUN = mean))
Figure4<-ggplot(DiaoCount,
                aes(x=movement, y=..density..))+
  geom_histogram(alpha=.7, binwidth = 1)+
  facet_wrap(~genotype, ncol=1)+
  coord_cartesian(xlim=c(0,40))+
  xlab("Number of cells GFP moves")+
  ylab("Density")+theme_bw()+
  geom_rug(aes(y=0),position = position_jitter(height = 0), alpha=0.4)+
  geom_vline(data = vline_df, aes(xintercept = Medians), linetype = 2, colour = "red4")+
  geom_vline(data = vline2_df, aes(xintercept = Means), linetype = 1, colour = "red4")+
  theme(text=element_text(size=25))


data1 <- DiaoCount[DiaoCount$genotype=="Col-0",]$movement
data2 <- DiaoCount[DiaoCount$genotype=="atfh2-1",]$movement
Figure5 <- boot_median(data1,data2)
Figure6 <- boot_mean(data1,data2)

DiaoLayer$genotype<-factor(DiaoLayer$genotype, levels(as.factor(DiaoLayer$genotype))[c(2,1)])
vline_df <- data.frame(genotype = levels(DiaoLayer$genotype),
                       Medians = tapply(X = DiaoLayer$movement, INDEX = DiaoLayer$genotype,
                                        FUN = median))
vline2_df <- data.frame(genotype = levels(DiaoLayer$genotype),
                        Means = tapply(X = DiaoLayer$movement, INDEX = DiaoLayer$genotype,
                                       FUN = mean))
Figure7<-ggplot(DiaoLayer,
                aes(x=movement, y=..density..))+
  geom_histogram(alpha=.7, binwidth = 1)+
  facet_wrap(~genotype, ncol=1)+
  coord_cartesian(xlim=c(0,5))+
  xlab("Number of layers GFP moves")+
  ylab("Density")+theme_bw()+
  geom_rug(aes(y=0),position = position_jitter(height = 0), alpha=0.4)+
  geom_vline(data = vline_df, aes(xintercept = Medians), linetype = 2, colour = "red4")+
  geom_vline(data = vline2_df, aes(xintercept = Means), linetype = 1, colour = "red4")+
  theme(text=element_text(size=25))


data1 <- DiaoLayer[DiaoLayer$genotype=="Col-0",]$movement
data2 <- DiaoLayer[DiaoLayer$genotype=="atfh2-1",]$movement
Figure8 <- boot_median(data1,data2)
Figure9 <- boot_mean(data1,data2)

lay <- rbind(c(-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0),
             c(1,1,1,1,3,3,3,3,5,5,5,5),
             c(2,2,2,2,4,4,4,4,6,6,6,6))

FigureA <-arrangeGrob(textGrob(""),left=textGrob("A)", x= unit(1, "npc"), y=unit(.95, "npc"),gp=gpar(fontsize=30)))
FigureB <-arrangeGrob(textGrob(""),left=textGrob("", x= unit(1, "npc"), y=unit(.95, "npc"),gp=gpar(fontsize=30)))
Figure1B<-arrangeGrob(Figure1, left=textGrob("B)", x= unit(1, "npc"), y=unit(.95, "npc"),gp=gpar(fontsize=30)))
Figure4B<-arrangeGrob(Figure4, left=textGrob("C)", x= unit(1, "npc"), y=unit(.95, "npc"),gp=gpar(fontsize=30)))
Figure7B<-arrangeGrob(Figure7, left=textGrob("D)", x= unit(1, "npc"), y=unit(.95, "npc"),gp=gpar(fontsize=30)))

# svglite("letterfig2.2.svg", width=20, height = 14)
grid.arrange(FigureA,
             FigureB,
             Figure1B,
             Figure2,
             Figure4B,
             Figure5,
             Figure7B,
             Figure9,
             layout_matrix = lay)
# dev.off()
