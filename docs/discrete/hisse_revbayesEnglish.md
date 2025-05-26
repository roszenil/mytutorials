---
title: HiSSE in Revbayes
layout: home
nav_order: 8
index: true
redirect: false
parent: Tutorials
math: katex
---
Created by Rosana Zenil-Ferguson for MOLE Workshop at Woods Hole, MA(May 2025). 
Data from di Stilio and Zenil-Ferguson. All traits considered. 2025. *Submitted*. 

For this tutorial we will create hidden-dependent state model for a pollinization with two states: Insect pollination (I) is denoted by 0 and wind pollination (W) is 1.

## Input

103 taxa in a bifurcating and ultrametric tree

 + List with species and their pollinization [here](files/poliniza_datos.csv)
 
 + Phylogenetic tree- [here](files/poliniza_arbol.tre)
 
 + The whole [HiSSE Rev code](files/hisse.Rev)
 
 ## RevBayes code
 
 1. Let's start with the vector ``moves`` that saves all the proposal moves for each of the parameters. Also the vector ``monitors`` saves the inference. Especially the posterior distribution for each of the parameters, so we can't forget that.
 
``` 
# Proposal moves (moves), and monitors for inference (the results of an MCMC)
moves = VectorMoves()
monitors = VectorMonitors()
```

2. The number of states is now different

```
NUM_STATES = 2
NUM_HIDDEN = 2
NUM_RATES = NUM_STATES * NUM_HIDDEN
```
 
3. Reading our data and the phylogeny

``` 
### Phylogenetic tree
observed_phylogeny <- readTrees("poliniza_arbol.tre")[1]


## Data
## 0 = Insect pollinated
## 1 = Wind pollinated
data <- readCharacterDataDelimited("poliniza_datos.csv",
stateLabels=2,
type="NaturalNumbers",
delimiter=",",
header=TRUE)

taxa <- observed_phylogeny.taxa()
```

4. But now we have to extend our 0 and 1 to 0A, 0B, 1A and 1B
```
data_exp <- data.expandCharacters( NUM_HIDDEN )
```

5. Some basic statistics of our data and tree

```
taxa <- observed_phylogeny.taxa()
root_age <- observed_phylogeny.rootAge()
```

6.  Prior distributions of the parameters

We will use the Gamma distribution as  **prior distribution** for the transition rates 

```
shape_pr := 0.5
rate_pr = observed_phylogeny.treeLength()/5
#### Create a Q-matrix size 4x4 
for (i in 1:NUM_RATES) {
for (j in 1:NUM_RATES) {
q[i][j] := 0.0
}
}

# In RevBayes 1=0A, 2=1A, 3=0B, and 4=1B
q_01A ~ dnGamma(shape=shape_pr, rate=rate_pr)
moves.append(mvScale(q_01A, weight=2 ))
q_10A ~ dnGamma(shape=shape_pr, rate=rate_pr)
moves.append(mvScale(q_10A, weight=2 ))
q_01B ~ dnGamma(shape=shape_pr, rate=rate_pr)
moves.append(mvScale(q_01B, weight=2 ))
q_10B ~ dnGamma(shape=shape_pr, rate=rate_pr)
moves.append(mvScale(q_10B, weight=2 ))
q[1][2] :=q_01A
q[3][4] :=q_01B
q[2][1] :=q_10A
q[3][4] :=q_10B
```

7. Now we need the transitions between the hidden states
```
hidden_rate1 ~ dnExponential(rate_pr)
moves.append(mvScale(hidden_rate1,lambda=0.2,tune=true,weight=2))
hidden_rate2 ~ dnExponential(rate_pr)
moves.append(mvScale(hidden_rate2,lambda=0.2,tune=true,weight=2))
#### Now drop them in the Q matrix
q[1][3] := hidden_rate1
q[2][4] := hidden_rate1
q[3][1] := hidden_rate2
q[4][2] := hidden_rate2
### Here we could define all of them different. This is the most traditional way. 
```

8. Build the CTMC called in comparative methods the Mk2 (Markov model with 2 states) using its essential tool: the Q-matrix

```
# The Q-matrix is an infinitesimal matrix, meaning it is a derivative of the probability matrix. 

rate_matrix := fnFreeK(q, rescaled=false, matrixExponentialMethod="scalingAndSquaring"))
```

9. Defining the diversification rates 

```
### Create prior parameters of the diversification rates
H = 0.5
rate_mean <- ln(ln(103/2.0) /root_age) # This is from Magallon and Sanderson (2001)
rate_sd <- 2*H

#We will be careful for the diversification rates of 0A y 0B. We will use a multiplier called sepciation_alpha for states in A and speciation_alpha+speciation_beta for B. We build them this way because A and B are hidden states derivated from the same main state so in theory they should be not so different. However, you can decide to define them separatedly. 

for (i in 1:NUM_STATES) {
### Create a lognormal distributed variable for the speciation rate
speciation_alpha[i] ~ dnNormal(mean=rate_mean,sd=rate_sd)
moves.append(mvSlide(speciation_alpha[i],delta=0.20,tune=true,weight=2.0))

### Create a lognormal distributed variable for the extinction rate
extinction_alpha[i] ~ dnNormal(mean=rate_mean,sd=rate_sd)
moves.append(mvSlide(extinction_alpha[i],delta=0.20,tune=true,weight=2.0))
}


for (i in 1:NUM_HIDDEN) {

### Create an exponential distributed variable for the speciation rate
speciation_beta[i] ~ dnExp(1.0)
moves.append(mvScale(speciation_beta[i],lambda=0.20,tune=true,weight=2.0))

### Create an normal distributed variable for the extinction
extinction_beta[i] ~ dnNormal(0.0,1.0)
moves.append(mvSlide(extinction_beta[i],delta=0.20,tune=true,weight=2.0))

}

for (j in 1:NUM_HIDDEN) {
for (i in 1:NUM_STATES) {
if ( j == 1) {
speciation[i] := exp( speciation_alpha[i] )
extinction[i] := exp( extinction_alpha[i] )
} else {
index = i+(j*NUM_STATES)-NUM_STATES
speciation[index] := exp( speciation_alpha[i] + speciation_beta[j-1] )
extinction[index] := exp( extinction_alpha[i] + extinction_beta[j-1] )
}
}
}

```

10. Root estimation

We do not know if the most common recent ancestor of all taxa were wind or insect pollinated, so we need to estimate the value of the root. In the Bayesian statistics framework the frequencies at the root are two extra parametes that we need to estimate. Therefore, we will assume initial equal frequencies for each state. Since there are two states we will need a bivariate prior distribution as the prior for a vector with two frequencies that add up to 1. A very useful distribution for this goal is the Dirichlet, a multivariate proability distribution that allows us to assign equal frequency to both states. 

``` 
root_frequencies ~ dnDirichlet(rep(1,NUM_STATES))

# Two proposals to explore those values at the root multiple times. 

moves.append(mvBetaSimplex(root_frequencies, alpha=0.5, weight=2))

moves.append(mvElementSwapSimplex(root_frequencies, weight=3))
```

11. Sampling fraction.
For diversification models we need to know how many lineages we have sampled out of the total lineages that are in the phylogenetic tree.
```
rho <- 103/200
```

12. Last step to build the HiSSE model

In order to connect all of the nodes in this model we have to merge them under a probability distribution for phylogenetic trees. This probability distribution is the one that internally calculates the likelihood function throughout the tree and considers the prior distribution for the parameters throughout the structure of the tree. 

```
#  As it name indicates the dnCDBDP is a  probability distribution creates a  Markov model that is Character Dependent Birth and Death Process over the phylogenetic tree

hisse ~ dnCDBDP( rootAge = root_age,
speciationRates   = speciation,
extinctionRates   = extinction,
Q                 = rate_matrix,
pi                = root_frequencies,
rho               = rho,
delta             = 1,
condition         = "time")
```



13. Calculating the likelihood function

Up to this point we have ignored the data. Of course, we need the data to calculate the likelihood. In RevBayes we do this throughout a function called ``clamp()``. This makes our  ``dnCDBDP()`` consider the values on the tips of the tree

``` 
# Clamp or data and tree to the model
hisse.clamp( observed_phylogeny )
hisse.clampCharData( data_exp ) #note the clamping on the expanded dataset

```

Once the data are clamped we will have a  **posterior distribution for the HiSSE model** and we can perform our inference

## Storing our Bayesian inference 

We have calculated the posterior distribution but now we need to focus on how to obtain the inference. If we look closely, throughout the code we specified  ``moves`` that are the proposals of the MCMC for each of the parameters. 

1. A brief description of the proposals

+ ``mvScale( q_01, weight=2 )`` this proposal rescales the original value of $$q_{01}$$ twice per iteration.
+ ``mvSlide(speciation_alpha[i],delta=0.20,tune=true,weight=2.0)`` this proposal is a window that slides through values of the parameter and tunes the size of the window as the MCMC advances. 
+ ``mvBetaSimplex(root_frequencies, alpha=0.5, weight=2)``  this function proposes two values in the interval (0,1) that added result in 1. They represent the frequency at which we would find the state 0 or the state 1 at the root of the tree.
+ ``mvElementSwapSimplex(root_frequencies, weight=3)`` this function proposes to swap the frequencies. For example, if we had (0.4, 0.6), elementswap exchanges the position of those values (0.6,0.4).

Just as discussed in the lectures, each of this proposals are going to contribute to the posterior odds, and if those improve they may get accepted.

2. The "store in a box" step of RevBayes

This is an important step for RevBayes software. We want to "store" the whole graphical object to be able to manipulate it. The function `model()` allows us to do this. 

``` 
# Mymodel is like a box that takes care of the whole graphical model
mymodel = model(rate_matrix)
```

3. The monitors follow the inferential process

The ``monitors`` store all our inference, without them you won't be able to retrieve your posterior distributions. There are many types of monitors as you will see below. 

``` 
## This monitor saves the whole posterior distribution for each of the parameters
monitors.append(mnModel(filename="output/hisse_pollination.log", printgen=1))

## This monitor prints in screen one parameter so you know that is running and how far.
monitors.append(mnScreen(printgen=10, q_01A, q_10A, speciation, extinction))

## This monitor tracks what is going on with the ancestral state reconstruction at each node
monitors.append(mnJointConditionalAncestralState(tree=hisse cdbdp=timetree, type="NaturalNumbers", printgen=1000, withTips=true, withStartStates=false, filename="output/anc_states_hisse_pollination.log"))

## This monitor creates stochastic character maps (evolution over the branches)- this is very slow for very complicated models

# We won't run the stochastic map here because it takes long time but this is how you would do it.
# monitors.append( mnStochasticCharacterMap(hisse,printgen=100,filename="output/stochmap_hisse_pollination.log", include_simmap=true))
```

+ ``mnModel`` Saves the samples of the posterior distribution created by the MCMC algorithm
+ ``mnScreen`` Prints in screen so you know something is happening
+ ``mnJointConditionalAncestralState`` Saves the ancestral state reconstruction utilizing the marginal posterior probability (be careful this is very different to what phytools or other software does) utilizando la probabilidad posterior marginal
+ ``mnStochasticCharacterMap`` calculates stochastic character maps that are transitions that occur along the branches. This is very important yet very difficult to calculate. 


4. Run your MCMC

**Remember to always run two chains to show that the MCMC has converge**
``` 
#### nruns=2, running this for 50,000 generations. 

mymcmc = mcmc(mymodel, monitors, moves, nruns=2, moveschedule="random")

### pre-burnin to tune the proposals 20% of the sample
#mymcmc.burnin(generations=2000,tuningInterval=100)

### run the MCMC
mymcmc.run(generations=150000)
```

5. Create summaries of your ancestral state reconstruction 

```
# Generating an marginal posterior distribution for each of the nodes to calculate their ancestral state reconstruction
###anc_state_trace = readAncestralStateTrace("output/anc_states_hisse_pollination_run_1.log")
###ancestralStateTree(tree=observed_phylogeny, ancestral_state_trace_vector=anc_state_trace, include_start_states=false, file="output/asr_hisse_polinizador.tree", summary_statistic="MAP", reconstruction="marginal")

q()
```

![](images/hissegraphical.png)
