# Bayesian Adaptive Direct Search

## Algorithm

BADS follows a mesh adaptive direct search (MADS) procedure that alternates **poll** steps and **search** steps (see **Fig 1**). 

- In the **poll** phase, points are evaluated on a (random) mesh by taking steps in one direction at a time, until an improvement is found or all directions have been tried. The step size is doubled in case of success, halved otherwise. 
- In the **search** step, a Gaussian process (GP) is fit to a (local) subset of the points evaluated so far. Points to evaluate during search are iteratively chosen by maximizing the *Expected Improvement* w.r.t. the current optimum [2,3].

**Fig 1: BADS procedure** ![BADS cartoon](https://github.com/lacerbi/bads/blob/master/figures/bads-cartoon.png "Fig 1: BADS cartoon")

Adherence to the MADS framework guarrantees convergence to a (local) stationary point of a noiseless function under general conditions [1]. The basic scheme is enhanced with heuristics to accelerate the poll step, to update the GP hyper-parameters, to generate a good set of candidate points in the search step, and to deal robustly with noisy functions.

## Methods

We tested the performance of BADS and of many optimization algorithms in MATLAB (**Fig 2**), such as Nelder-Mead (`fminsearch`), active-set SQP (`fmincon`), pattern search (`patternsearch`), Multi-Level Coordinate Search (`mcs`), simulated annealing (`simulannealbnd`), genetic algorithms (`ga`), particle swarm (`particleswarm`), CMA-ES (`cmaes`). 
In particular, MCS and CMA-ES are state-of-the-art methods for non-convex, derivative-free optimization [4]. For all algorithms, including BADS, we used default settings (no fine-tuning).

First, we benchmarked on standard test functions (**Fig 2**):
- noiseless (*cec14* test suite; 8 functions defined on `D = 2,3,5,10,20`); and
- noisy (*cec14-noisy*; as before plus random Gaussian noise at each evaluation). 

Then, we benchmarked on actual model-fitting problems from four projects in our lab (**Fig 3**), where the objective to optimize is the log-likelihood as a function of model parameters, for a given model and dataset (noiseless: *loc*, *causinf*; noisy: *change*, *wrm*; 4-6 distinct test functions per problem).

For the benchmark, we ran 50 independent runs of each algorithm on each test function, with randomized starting points and a budget of `500 * D` function evaluations. If an algorithm terminated before depleting the budget, it was restarted from a new random point.
We consider a run *successful* if the best achieved function value is within *eps* from the true optimum `f_opt` (or our best estimate thereof). We set a target of *eps* = 0.1 for noiseless functions and *eps* = 1 for noisy functions (smaller differences in log-likelihood are negligible due to data and model uncertainty).

## Results



## Conclusions

BADS combines MADS [1] and Bayesian Optimization [2], similarly in principle to [3], but in a novel approach that makes it effective for model fitting. Previous work in Bayesian Optimization targets only slow simulations on much longer timescales [2,3]. 
BADS can be used out-of-the-box and will be made freely available as a MATLAB package. These features make BADS of outstanding interest for modellers in computational neuroscience.

## References

1. Audet C & Dennis Jr JE (2006). "Mesh Adaptive Direct Search Algorithms for Constrained Optimization", *SIAM J Opt* 
([link](http://epubs.siam.org/doi/abs/10.1137/040603371)).
2. Snoek *et al*. (2012). "Practical Bayesian Optimization of Machine Learning Algorithms", *NIPS* 2012(4522) ([link](https://papers.nips.cc/paper/4522-practical-bayesian-optimization-of-machine-learning-algorithms)).
3. Gramacy RB & Le Digabel S (2015). "The mesh adaptive direct search algorithm with treed Gaussian process surrogates", *Pac J Opt* ([link](http://www.optimization-online.org/DB_HTML/2011/07/3090.html)).
4. Rios LM & Sahinidis NV (2013). "Derivative-free optimization: a review of algorithms and comparison of software implementations", *J Global Opt* ([link](http://link.springer.com/article/10.1007/s10898-012-9951-y)).

