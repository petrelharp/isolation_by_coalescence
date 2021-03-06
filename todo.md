
1. Rewrite SLiM sims to have larger $\sigma_d$ and bigger genomes.
2. Run commute time inference on SLiM results.
3. Run commute time inference on tortoises.
4. Read and cite Hanks and Hooten.
5. Note that circuit theory assumes symmetry, ensures identifiability
6. Note that resistance and coalescence are the same with isotropic graphs (mcrae2006)
7. Make point that we're no longer data-limited, unlike the reams of papers using IBR on 10 microsats.

# Literature notes

hanks2017modeling -- discusses asymmetry, with application to trout in a stream network
    Has a mixed model where the driving Gaussian noise has the covariance function 
    of the stationary distribution of a discrete-space-time random field.
    The random field is a discrete model of abundance with sources, sinks, and migration -
    reaction-diffusion with diffusion from a MC and reaction equal to births minus deaths;
    and modeling the reaction term as spatial white noise.
    No reasons are given why the covariance function for space-time abundance should equal
    the covariance function underlying allelic identity.


## Resistance

McRae 2006 - models linearized Fst (which is Fst / (1-Fst) = (pi[i,j] - pi[i,i])/2*pi[i,i]);
    cites Slatkin to say it's a difference of coalescence times;
    and then on an isotropic graph coalescence is one-quarter commute times.
    Says that in nonisotropic cases it's a convenient approximation
    that is more computationally efficient.
    Computes coalescence times and compares these to resistance distance
    for a number of cases (all symmetric);
    shows correlation is not perfect.
    Compares LCP and IBD to IBR.

Cushman 2006 (cushman2006complex) ("Gene Flow in Complex Landscapes: Testing Multiple Hypotheses With
    Causal Modeling") - "isolation by resistance hypothesis"
    applied to black bears: uses mean divergence with 9 microsats
    and compares 108 different landscape resistance hypotheses
    (factorial combinations of four layers)

Hanks and Hooten 2013 (hanks2013circuit) - Does Bayesian inference with resistance distance,
    connects to Intrinsic Conditional AutoRegressive models (ICAR; Besag 1974).
    Has another method to compute R.
    From covariance matrix S, defines distance D_ij = Sii + Sjj - 2 Sij;
    and then defines a GMRF whose covariance matrix induces D,
    with noise modeled as Wishart.

Graves et al 2013 (graves2013current) "Current approaches using genetic distances produce poor estimates of landscape resistance to interindividual dispersal"
    Pyrrhic victory: gets best Mantel-correlated resistance parameters,
    but the Mantel surface is very flat. Used "Dps" for genetic distance, 
    but cites a bunch of other papers saying that others are highly correlated.

## Other migration inference methods

Wilson and Rannala 2003 (wilson2003bayesian) ("Bayesian Inference of Recent Migration Rates Using Multilocus Genotypes,")
    This is BayesAss; does a coalescent theory model, Bayesian.

Some of the MIGRATE papers:
    and Beerli and Felsenstein 1999 (beerli1999maximum) -- also coalescent theory; does only two populations -- does better than using Fst = 1/4nm
    Beerli 2006 ("Comparison of Bayesian and Maximum-Likelihood Inference of Population Genetic Parameters")
    Beerli and Palczewski 2010 ("Unified Framework to Evaluate Panmixia and Migration Direction Among Multiple Sampling Locations") -- thermodynamic integration

Hey and Nielsen 2007 (hey2007integration) -- this is IMa;
    does a structured coalescent theory model.

Hey and Nielsen 2004 (hey2004multilocus) -- this is one of the IM papers; structured coalescent model.

## Examples

(see Hanks and Hooten 2013 for a bunch)

Wang, Savage, and Shaffer (wang2009landscape) - least cost paths, tiger salamander

McRae and Beier 2007 (mcrae2007circuit,) ("Circuit Theory Predicts Gene Flow in Plant
    and Animal Populations") - compares LCP and IBD to IBR in two empirical cases.


Rioux Paquette et al 2014 (riouxpaquette2014modelling) ("Modelling the Dispersal of the Two Main Hosts of the Raccoon Rabies")
    An example, again using microsats.
    Comes to some strange conclusions, like different resistance surfaces for males and females,
    based on between-male and between-female distances (!?!?).
