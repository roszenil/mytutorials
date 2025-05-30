---
title: "HiSSE visualization"
layout: home
nav_order: 9
index: true
redirect: false
parent: Tutorials
math: katex
---

Created by Rosana Zenil-Ferguson for MOLE workshop, Woods Hole, MA (May 2025)

## Setting up in your computer

There are two options to run smoothly this tutorial

### Option 1

1. Make sure you have the most up-to-date R and RStudio in your personal computer.
2. Make sure you create a project by doing Files> Open new project > New directory > New project
3. Name your project a something memorable like discrete_mole2025 and save it in a location where it will be safe. See image below
![](images/rproject.png)
4. Download all the files in this tutorial in that exact folder and everything should work smoothly

### Option 2

1. Open this [link](https://posit.cloud/content/10437789) 
2. Create a free PositCloud account
3. You have the full project!

## Library organization

We need the following packages for this tutorial.  **Important**: We need an specific version of RevGadgets and not the version currently in CRAN to do the stochastic mapping. 

```
# install.packages("devtools")
# install.packages("coda")
# install.packages("phytools")
# install.packages ("ggplot2")
# install.packages ("dplyr")
## Instalar treeio
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##ggtree
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# Do not do this: install.packages("RevGadgets") 
# We want the version below devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)


library(coda)
library(phytools)
library(ggplot2)
library(dplyr)
library(RevGadgets)
```

## MCMC convergence

Check for convergence with the coda library

```
# Make sure you are where your files are supposed to be
setwd()

# With coda
mcmc_run1 <- readTrace(path ="hisse_pollination_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)

# If you run into issues use it as below
# mcmc_trace[[1]] <- as.mcmc(mcmc_trace[[1]]
```

Better job with convergence plots

```
# ggplot+revgadgets
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

# Make sure you have cut sufficiently the burn-in
```


## Rate estimates

We are going to plot the posterior distribution of the parameters. For the HiSSE model we have a lot of them 

+ Transition rates between main states $$q_{01A}, q_{10A}, q_{01B}, q_{10B}$$
+ Transition rates between hidden states $$\alpha, \beta$$
+ Net diversification rates $$r_{0A}=\lambda_{0A}-\mu_{0A},r_{1A}=\lambda_{1A}-\mu_{1A},r_{0B}=\lambda_{0B}-\mu_{0B}, r_{1B}=\lambda_{1B}-\mu_{1B}$$
+ Alternatively: Extinction fractions, turnover $$\epsilon_{0A}=\mu_{0A}/\lambda_{0A}, \epsilon_{1A}=\mu_{1A}/\lambda_{1A}, \epsilon_{0B}=\mu_{0B}/\lambda_{0B},\epsilon_{1B}=\mu_{1B}/\lambda_{1B}$$
+ Root frequencies

For more complex models like HiSSE I always prefer to plot them with violins to have more control over colors and overlapping. So I do this with ggplot2 but you can do posterior densities with RevGadgets.

## Transition rate plots

```
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


## Hidden state transitions

hidden_rates<- data.frame(dens=c(hisse$hidden_rate1,hisse$hidden_rate2) ,rate=rep(c("alpha","beta"),each=length(hisse$hidden_rate1)))

violin_hidden<- ggplot(hidden_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Hidden state transitions ")+
  scale_fill_manual( values = traitcols[1:2])+
  xlab("Transition type")+
  ylab("Rate")+
  theme_classic()
violin_hidden
```

## Net diversification rates


```
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
```

**What can we quickly conclude from this plot?**

## Hypothesis testing: Is pollinization linked to diversification?

We build the test statistics for the differences between net diversifications. We will compare differences between 0 and 1 for each hidden state
+ $$T_A= (\lambda_{0A}=\mu_{0A})-(\lambda_{1A}=\mu_{1A})$$ 
+ $$T_B= (\lambda_{0B}=\mu_{0B})-(\lambda_{1B}=\mu_{1B})$$ 

If these differences have probability larger than 5% of being zero then this means diversification is not correlated to our traits. 

```
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
# 0 crosses these posterior distributions, but we have to check the probability.

```

Formalizing the hypothesis test
$$H_0: T_A=0$$ and  $$T_B=0$$

```

# dplyr

quantile_diff<- T_diff %>% group_by(difference)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
quantile_diff
```

Let's observe that the credibility interval at 95% for $$T_A$$ is (-0.219, 0.425) and for $$T_B$$ is (-0.554, 0.0766). Since 0 belongs to these two intervals then $$P(H_0\lvert Datos)>0.05$$, this means that the tempo of diversification for 0 and 1 is equal. **Important**: Note that in a Bayesian framework we do not use p-values, we don't use the word *significance* or *rejection*. We don't use those words because we are using probability distributions and not likelihood. Be careful when writing your results.

**Second part**: What about the hidden states?

We build test statistics for the net diversification differences between A and B for state 0 and 1 respectively. 

+ $$T_0= (\lambda_{0A}=\mu_{0A})-(\lambda_{0B}=\mu_{0B})$$ 
+ $$T_1= (\lambda_{1A}=\mu_{1A})-(\lambda_{1B}=\mu_{1B})$$ 

If these differences are 0 with probability larger than 0.05 then we have that there are changes in the tempo of diversification due to something else that we didn't directly measure. 

```
#ggplot2
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

# Observemos que 0 cruza la diferencia entre las tasas lo que quiere decir que pueden ser iguales
```

Formalizing the hypothesis testing
$$ H_0: T_0=0$$ and  $$T_1=0$$

```
# dplyr

quantile_diff <- T_diff %>% group_by(difference)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
quantile_diff
```

We see that the credible interval at 95% for $$T_0$$ is (-0.488, -0.214) and for $$T_1$$ es (-0.989, -0.321). Then, zero **does not** belong to these intervals then $$P(H_0\lvert Datos)< 0.05$$. This means, that the tempo of diversification is correlated to the hidden states. This result along with the previous one indicate that the correct model of state-dependent diversification for this system is the CID-2. 


### Plot ancestral states map 

+ Dowload the marginal ancestral state reconstruction from RevBayes [here](files/asr_hisse_polinizador)

Let's look at what happened

```
anc_states <- processAncStates(path ="asr_hisse_polinizador.tree",state_labels=c("0"="insect A","1"="wind A","2"="insect B", "3"="wind B"))
plotAncStatesMAP(t = anc_states, tree_layout="rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2,
                 node_color = c("#D9081D","#F8D3D8","#072AC8","#A2D6F9"))
```

