
# install.packages("devtools")
# install.packages("coda")
# install.packages("phytools")
# install.packages ("ggplot2")
# install.packages ("dplyr")
## Instalar treeio
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##Instalar ggtree
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# Importante: no instales install.packages("RevGadgets") porque esto instala desde CRAN. Queremos una version especifica. Por favor sigue la instruccion de abajo. Si ya habias instalado revgadgets por favor reinicia y reinstala asi

# devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)


library(coda)
library(phytools)
library(ggplot2)
library(dplyr)
library(RevGadgets)


# Make sure you are in the correct working directory or in the correct project
#setwd( )

# Package code
mcmc_run1 <- readTrace(path ="bisse_pollination_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)


Better job with convergence plot

#  ggplot2 y con RevGadgets
mcmc_run1 <- readTrace(path ="~/Dropbox/Teaching/Workshops/UNAM2025/discrete_trait/files/bisse_pollination_run_1.log", burnin = 0.1)
mcmc_run1<-data.frame(mcmc_run1)
mcmc_run1<- cbind(mcmc_run1,run=rep("run 1",length(mcmc_run1$Iteration)))
mcmc_run2 <- readTrace(path ="~/Dropbox/Teaching/Workshops/UNAM2025/discrete_trait/files/bisse_pollination_run_2.log", burnin = 0.1)
mcmc_run2<-data.frame(mcmc_run2)
mcmc_run2<- cbind(mcmc_run2,run=rep("run 2",length(mcmc_run2$Iteration)))

mcmc_table<-rbind(mcmc_run1,mcmc_run2)

trace_plot<- ggplot(mcmc_table, aes(x=Iteration,y=Posterior,group=run))+
  geom_line(aes(color=run))+
  theme_classic()
trace_plot



traitcols<-c("#3D348B","#F18701")

bisse<- read.table("bisse_pollination_run_1.log", header=TRUE)
bisse<- bisse[-seq(1,5000,1),] # burn in removal

transition_rates<- data.frame(dens=c(bisse$q_01, bisse$q_10) ,rate=rep(c("q_01","q_10"),each=length(bisse$q_01)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Transition rates")+
  scale_fill_manual( values = traitcols)+
  xlab("Transition type")+
  ylab("Rate")+
  theme_classic()
violin_transitions


divcols<-c("#E63946","#1D3557")
netdiversification_rates<- data.frame(dens=c(bisse$net_diversification.1., bisse$net_diversification.2.) ,rate=rep(c("r_0","r_1"),each=length(bisse$net_diversification.1.)))

violin_diversification<- ggplot(netdiversification_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Net diversification rates")+
  scale_fill_manual( values = divcols)+
  xlab("")+
  ylab("Rate")+
  theme_classic()
violin_diversification


freqcols<-c("hotpink","darkblue")
rootfreq<- data.frame(dens=c(bisse$root_frequencies.1., bisse$root_frequencies.2.) ,rate=rep(c("state 0 freq","state 1 freq"),each=length(bisse$root_frequencies.1.)))

violin_rootfreq<- ggplot(rootfreq,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Root frequencies")+
  scale_fill_manual( values = freqcols)+
  theme_classic()
violin_rootfreq



difcols<-c("#1DB32C","#BFEEC3")
T_diff<- data.frame(dens=(bisse$net_diversification.1.-bisse$net_diversification.2.),difference=rep("T",each=length(bisse$net_diversification.1.)))

violin_difference<- ggplot(T_diff,aes(x=difference,y=dens, fill=difference))+
  geom_violin(trim=FALSE)+
  labs(title="Summary statistic")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

# Let's not that zero crosses this violin in the tail of the distribution. However, we want to know if zero belongs to the credible interval
violin_difference

## Does zero belong to the credible interval?
quantile ((bisse$net_diversification.1.-bisse$net_diversification.2.),probs=c(0.025,0.975))
#Answer: NO!
