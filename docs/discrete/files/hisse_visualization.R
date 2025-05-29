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


# Check that you are in the correct project or wd
#setwd()

# CODA
mcmc_run1 <- readTrace(path ="hisse_pollination_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)

#  ggplot2 and RevGadgets
mcmc_run1 <- readTrace(path ="hisse_pollination_run_1.log", burnin = 0.1)
mcmc_run1<-data.frame(mcmc_run1)
mcmc_run1<- cbind(mcmc_run1,run=rep("run 1",length(mcmc_run1$Iteration)))

mcmc_run2 <- readTrace(path ="hisse_pollination_run_2.log", burnin = 0.1)
mcmc_run2<-data.frame(mcmc_run2)
mcmc_run2<- cbind(mcmc_run2,run=rep("run 2",length(mcmc_run2$Iteration)))

mcmc_table<-rbind(mcmc_run1,mcmc_run2)

trace_plot<- ggplot(mcmc_table, aes(x=Iteration,y=Posterior,group=run))+
  geom_line(aes(color=run))+
  theme_classic()
trace_plot


## Transition rates
# ggplot2

traitcols<-c("#3D348B","#7678ED","#F18701", "#F35B04")

hisse<- read.table("hisse_pollination_run_1.log", header=TRUE)
hisse<- hisse[-seq(1,15000,1),] # make sure you are cutting the burn in!

transition_rates<- data.frame(dens=c(hisse$q_01A,hisse$q_01B, hisse$q_10A, hisse$q_10B) ,rate=rep(c("q_01A","q_01B","q_10A","q_10B"),each=length(hisse$q_01A)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Transition Rates")+
  scale_fill_manual( values = traitcols)+
  xlab("Transition type")+
  ylab("Rate")+
  theme_classic()
violin_transitions

## Hidden rate transitions

hidden_rates<- data.frame(dens=c(hisse$hidden_rate1,hisse$hidden_rate2) ,rate=rep(c("alpha","beta"),each=length(hisse$hidden_rate1)))

violin_hidden<- ggplot(hidden_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Hidden state transitions ")+
  scale_fill_manual( values = traitcols[1:2])+
  xlab("Transition type")+
  ylab("Rate")+
  theme_classic()
violin_hidden



## Diversification rates

divcols<-c("#E63946","#F3A5AB","#1D3557", "#457B9D")

# In RevBayes 1=0A, 2=1A, 3=0B, and 4=1B
netdiversification_rates<- data.frame(dens=c(hisse$speciation.1.-hisse$extinction.1.,hisse$speciation.2.-hisse$extinction.2., hisse$speciation.3.-hisse$extinction.3.,hisse$speciation.4.-hisse$extinction.4.) ,rate=rep(c("net_div_0A","net_div_1A","net_div_0B","net_div_1B"),each=length(hisse$speciation.1.)))

violin_diversification<- ggplot(netdiversification_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Tasas de transicion")+
  scale_fill_manual( values = divcols)+
  xlab("")+
  ylab("Rate")+
  theme_classic()
violin_diversification

  
  ## Hypothesis testing for diversification correlated to states
  
  # $$T_A= (\lambda_{0A}=\mu_{0A})-(\lambda_{1A}=\mu_{1A})$$ 
  # $$T_B= (\lambda_{0B}=\mu_{0B})-(\lambda_{1B}=\mu_{1B})$$ 
  

# ggplot2

difcols<-c("#FF006E","#FFC2DC")
T_diff<- data.frame(dens=c((hisse$speciation.1.-hisse$extinction.1.)-(hisse$speciation.2.-hisse$extinction.2.),(hisse$speciation.3.-hisse$extinction.3.)-(hisse$speciation.4.-hisse$extinction.4.)),difference=rep(c("T_A","T_B"),each=length(hisse$speciation.1.)))

violin_difference<- ggplot(T_diff,aes(x=difference,y=dens, fill=difference))+
  geom_violin(trim=FALSE)+
  labs(title="Test statistics")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  xlab("Test statistic")+
  ylab("Difference")+
  theme_classic()

violin_difference



# dplyr

quantile_diff<- T_diff %>% group_by(difference)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
quantile_diff

## Hypothesis testing for the amount of heterogeneity

#  + $$T_0= (\lambda_{0A}=\mu_{0A})-(\lambda_{0B}=\mu_{0B})$$ 
#  + $$T_1= (\lambda_{1A}=\mu_{1A})-(\lambda_{1B}=\mu_{1B})$$ 
  

# ggplot2

difcols<-c("#1DB32C","#BFEEC3")
T_diff<- data.frame(dens=c((hisse$speciation.1.-hisse$extinction.1.)-(hisse$speciation.3.-hisse$extinction.3.),(hisse$speciation.2.-hisse$extinction.2.)-(hisse$speciation.4.-hisse$extinction.4.)),difference=rep(c("T_0","T_1"),each=length(hisse$speciation.1.)))

violin_difference<- ggplot(T_diff,aes(x=difference,y=dens, fill=difference))+
  geom_violin(trim=FALSE)+
  labs(title="Test statistics")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  xlab("Test statistic")+
  ylab("Differences")+
  theme_classic()

violin_difference


quantile_diff <- T_diff %>% group_by(difference)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
quantile_diff
