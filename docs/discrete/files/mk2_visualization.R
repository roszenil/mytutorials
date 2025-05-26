
install.packages("devtools")
install.packages("coda")
install.packages("phytools")
install.packages ("ggplot2")

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("treeio")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)
library(coda)
library(RevGadgets)
library(phytools)
library(ggplot2)

mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)


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


mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
summarizeTrace(trace = mcmc_run1, vars =c("q_01","q_10","root_frequencies[1]","root_frequencies[2]"))

plotTrace(trace = mcmc_run1, vars = c("q_01","q_10"))[[1]]

plotTrace(trace = mcmc_run1, vars = c("root_frequencies[1]","root_frequencies[2]"))[[1]]



traitcols<-c("#7678ED","#F35B04")

mk2 <- read.table("mk2_polinizador_run_1.log", header=TRUE)
mk2<- mk2[-seq(1,5000,1),] # burn in

transition_rates<- data.frame(dens=c(mk2$q_01,mk2$q_10),rate=rep(c("q_01","q_10"),each=length(mk2$q_01)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Transition rates")+
  scale_fill_manual( values = traitcols)+
  theme_classic()
violin_transitions

D<- data.frame(dens=(mk2$q_01-mk2$q_10),rate=rep(c("Difference"),length(mk2$q_01)))

violin_difference<- ggplot(D,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Prueba de hipótesis")+
  scale_fill_manual( values = "hotpink")+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

violin_difference


diff_quantile <-quantile(D$dens, probs=c(0.025,0.975))
diff_quantile


anc_states <- processAncStates(path ="asr_mk2.tree",state_labels=c("0"="insect","1"="wind"))
plotAncStatesMAP(t = anc_states, tree_layout="rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2)

plotAncStatesPie(t=  anc_states,
                 state_transparency = 0.8,
                 node_labels_size=  5) 

mycolors= setNames(c("blue","darkorange"),c("0","1"))
file <- "poliniza_arbol.nex" #Has to be the nexus file as given by RevBayes
pol.tree <- readTrees(paths = file)[[1]][[1]]

mapsfile <- "stochmap_mk2_polinizador_run_1.log" 

stoch_map_df <- processStochMaps(pol.tree,
                                 mapsfile, 
                                 states = as.character(0:1), 
                                 burnin = 0.1)

plotStochMaps(tree=pol.tree,maps=stoch_map_df,tip_labels_size=0.5,colors=mycolors)
