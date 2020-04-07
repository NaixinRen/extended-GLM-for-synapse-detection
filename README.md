# Synapse Detection from Large-scale Spike Recording
 
Detecting synaptic connections using large-scale extracellular spike recordings presents a statistical challenge. Here we develop an extension of a Generalized Linear Model that describes the cross-correlograms between pairs of neurons to detect synaptic connections from large-scale spike recording. The extended GLM explicitly separates the cross-correlograms into fast synaptic effects and slow background fluctuations, and also incorporates two structural constraints learned from the whole network: presynaptic neuron type and the relationship between the synaptic latency and distance between pre- and postsynaptic neurons. 

[Ren N, Ito S, Hafizi H, Beggs JM, Stevenson IH (2020) Model-based detection of putative synaptic connections from spike recordings with latency and type constraints. *bioRxiv*](https://www.biorxiv.org/content/10.1101/2020.02.12.944496v1)

## Running the demo

demo.m shows an example of running the model on a simulated dataset.

### Core functions
* generate_correlogram() converts spike trains `spikes` into cross-correlograms `CCG`
```
[CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index); 
```
* learning_basis() learns a set of smooth basis functions `X` to capture the slow fluctuations in correlograms
```
X = learning_basis(CCG,ignore); 
```
* extendedGLM() fits the model on all the cross-correlogrms for each presynaptic neuron
```
model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));  
```
* detect_cnx() returns the final results based on the model fits
```
results = detect_cnx(model_fits,ignore,threshold);
```
where `threshold` is optional. 

### Results visualization

Compare overall connectivity matrices

![connectivity matrices](https://github.com/NaixinRen/extended-GLM-for-synapse-detection/raw/master/pics/ConnectivityMatrices_50neurons.png)

