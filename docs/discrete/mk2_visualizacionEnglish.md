---
title: "Results Mk2"
layout: home
nav_order: 4
index: true
redirect: false
parent: Temario
math: katex
---

Created by Rosana Zenil-Ferguson (May 2025)

## Organizing libraries

We need the following R packages. **Important**: The RevGadgets packages is available with ``install.packages("RevGadgets")``. However we need a special version for the stochastic maps that is not available on CRAN. If you don't want to do stochastic maps you can use the regular version. 

```
# install.packages("devtools")
# install.packages("coda")
# install.packages("phytools")
# install.packages ("ggplot2")

## Instalar treeio
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##Instalar ggtree
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# Important: Do not do install.packages("RevGadgets"). Do the version below

# devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)
library(coda)
library(RevGadgets)
library(phytools)
library(ggplot2)
```

## MCMC convergence

Check basic convergence

Download the following files
+ [First run](https://ixchelgzlzr.github.io/filo_bayes_UNAM/docs/discrete/files/mk2_polinizador_run_1.log)
+ [Second run(https://github.com/ixchelgzlzr/filo_bayes_UNAM/blob/main/docs/discrete/files/mk2_polinizador_run_2.log)

```
# Make sure you are in the correct working directory
#setwd("~/Teaching/Workshops/UNAM2025/discrete_trait/files")

# Coda package checks for basic convergence
mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)
```

Better job at visualizing convergence

```
# ggplot+ revGadgets
mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
mcmc_run1<-data.frame(mcmc_run1)
mcmc_run1<- cbind(mcmc_run1,run=rep("run 1",length(mcmc_run1$Iteration)))

mcmc_run2 <- readTrace(path ="mk2_polinizador_run_2.log", burnin = 0.1)
mcmc_run2<-data.frame(mcmc_run2)
mcmc_run2<- cbind(mcmc_run2,run=rep("run 2",length(mcmc_run2$Iteration)))

mcmc_table<-rbind(mcmc_run1,mcmc_run2)

trace_plot<- ggplot(mcmc_table, aes(x=Iteration,y=Posterior,group=run))+
              geom_line(aes(color=run))+
              theme_classic()
trace_plot

```

## Posterior distributions of rates

Remember that for Mk2 we estimated the following four parameters

+ $$q_{01}$$ - transition rate
+ $$q_{10}$$ - transition rate
+ Root frequencies


```
# revgadgets

mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
# Summary statistics
## MAP= maximum a posteriori
## Mean= mean of the posterior distribution
## Quantiles 2.5 and 97.5 define the credible interval
summarizeTrace(trace = mcmc_run1, vars =c("q_01","q_10","root_frequencies[1]","root_frequencies[2]"))

## Transition rates posterior distributions
plotTrace(trace = mcmc_run1, vars = c("q_01","q_10"))[[1]]

plotTrace(trace = mcmc_run1, vars = c("root_frequencies[1]","root_frequencies[2]"))[[1]]
```

Personally I prefer to graph violins, especially when there are more than two states.

```
# We do this with ggplot2

traitcols<-c("#7678ED","#F35B04") #Check out coolors.co for color choosing

mk2 <- read.table("mk2_polinizador_run_1.log", header=TRUE)
mk2<- mk2[-seq(1,5000,1),] # Never forget to cut your burn-in

transition_rates<- data.frame(dens=c(mk2$q_01,mk2$q_10),rate=rep(c("q_01","q_10"),each=length(mk2$q_01)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Transition rates")+
  scale_fill_manual( values = traitcols)+
  theme_classic()
violin_transitions

```

## Hypothesis testing: Are the transition rates between states equal?

One of the advantages of using posterior distribution is that we can quickly assess which parameters might be equal or different. However, the formal approach is to present a test statistic about the diference. How do we formalize this test in a Bayesian framework?

We build a summary statistic $$D= q_{01}-q_{10}$$ because difference and we plot it.

```
# Using ggplot2

D<- data.frame(dens=(mk2$q_01-mk2$q_10),rate=rep(c("Difference"),length(mk2$q_01)))

violin_difference<- ggplot(D,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Test statistic. Difference between rates")+
  scale_fill_manual( values = "hotpink")+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

# Zero is crossing the test statistic high probability values
violin_difference

```
Formalizing the hypothesis
$$ H_0: D=0$$


```
hpd <-quantile(D$dens, probs=c(0.025,0.975))
hpd
```

We observa that the highest posterior probability interval, known as the credible interval for $$D$$ with 95% of the probability of the posterior is (-0.19, 0.01). The value zero belongs to the credible interval. Therefore, the probability of the null hypothesis is larger than 5%, in mathematical terms $$P(H_0\lvert Data)>0.05$$. This means that the transition rates are equal with a probability higher than 5%.  **Important**: Note that in Bayesian hypothesis testing we don't use p-values, or significance, we do not reject or fail to reject, we simply state how **probable** a hypothesis is. 

## Ancestral reconstructions

There are two tools that help us reconstruct the past. These tools are:

1. Marginal ancestral state reconstruction. In a Bayesian framework, we estimate ancestral states by calculating the marginal probability of each node. This means that we calculate the probability of each state value at a single node, while "averaging" (integrating/summing) over the rest ot the nodes. 

    + [Phylogenetic tree with marginal ancestral state reconstruction](files/asr_mk2.tree)

    ```
       # We use RevGadgets

       anc_states <- processAncStates(path ="asr_mk2.tree",state_labels=c("0"="insect","1"="wind"))
       plotAncStatesMAP(t = anc_states, tree_layout="rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2)
       # produce the plot object, showing MAP states at nodes.
       # color corresponds to state, size to the state's posterior probability

    ```


2. Stochastic maps: These are simulations from root to tip of probable histories that happened throughout branches of the tree. Better than only reconstructing over nodes, but more difficult to calculate. We have multiples of this simulations but one way to summarize is to choose the maximum *a posterior* of all those simulations as we have it below. 

    + [Phylogenetic tree in Nexus format](files/poliniza_arbol.nex)
    + [Stochastic maps](files/stochmap_mk2_polinizador_run_1.log)

    ```
      #This is in development in RevGadgets hence the request to download the specific version above. 

      mycolors= setNames(c("blue","darkorange"),c("0","1"))

      # Read the tree
      file <- "poliniza_arbol.nex" #Has to be the nexus file as given by RevBayes
      pol.tree <- readTrees(paths = file)[[1]][[1]]

      # Read the stochastic maps
      mapsfile <- "stochmap_mk2_polinizador_run_1.log" 

      # Summarize them with the MAP
      stoch_map_df <- processStochMaps(pol.tree,
                                 mapsfile, 
                                 states = as.character(0:1), 
                                 burnin = 0.1)

      # Plot
      plotStochMaps(tree=pol.tree,maps=stoch_map_df,tip_labels_size=0.5,colors=mycolors)
    ```
