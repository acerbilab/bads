# Bayesian Adaptive Direct Search

Bayesian Adaptive Direct Search (BADS) is a novel Bayesian Optimization algorithm that achieves competitive performance at a small computational cost [1]. BADS has been intensively tested for fitting behavioral, cognitive, and neural models and is currently being used in more than a dozen projects in [our lab](http://www.cns.nyu.edu/malab/).

In our benchmark with artificial and real problems, only BADS performs consistently well on all problems: BADS generally outperforms all other methods on non-smooth or noisy functions with dimension `3 <= D <= 20`, and is competitive with methods such as `fmincon` (MATLAB) and Multilevel Coordinate Search on smooth, noiseless functions.

Code to run BADS out-of-the-box will be released here, together with a white paper (ETA: end of March 2017).

Please contact luigi.acerbi@nyu.edu for more information.


## Algorithm

BADS follows a mesh adaptive direct search (MADS) procedure [2] that alternates **poll** steps and **search** steps (see **Fig 1**). 

- In the **poll** phase, points are evaluated on a (random) mesh by taking steps in one direction at a time, until an improvement is found or all directions have been tried. The step size is doubled in case of success, halved otherwise. 
- In the **search** step, a Gaussian process (GP) is fit to a (local) subset of the points evaluated so far. Points to evaluate during search are iteratively chosen by maximizing the *Expected Improvement* w.r.t. the current optimum [3,4].

**Fig 1: BADS procedure** ![BADS procedure](https://github.com/lacerbi/bads-dev/blob/master/docs/figures/fig1-demo.png "Fig 1: BADS procedure")

Adherence to the MADS framework guarrantees convergence to a (local) stationary point of a noiseless function under general conditions [1]. The basic scheme is enhanced with heuristics to accelerate the poll step, to update the GP hyper-parameters, to generate a good set of candidate points in the search step, and to deal robustly with noisy functions.

## Methods

We tested the performance of BADS and of many optimization algorithms in MATLAB (**Fig 2**), such as Nelder-Mead (`fminsearch`), active-set SQP (`fmincon`), pattern search (`patternsearch`), Multi-Level Coordinate Search (`mcs`), simulated annealing (`simulannealbnd`), genetic algorithms (`ga`), particle swarm (`particleswarm`), CMA-ES (`cmaes`). 
In particular, MCS and CMA-ES are state-of-the-art methods for non-convex, derivative-free optimization [5]. For all algorithms, including BADS, we used default settings (no fine-tuning).

First, we benchmarked on standard test functions (**Fig 2**):
- noiseless (*cec14* test suite; 8 functions defined on `D = 2,3,5,10,20`); and
- noisy (*cec14-noisy*; as before plus random Gaussian noise at each evaluation). 

Then, we benchmarked on actual model-fitting problems from four projects in our lab (**Fig 3**), where the objective to optimize is the log-likelihood as a function of model parameters, for a given model and dataset (noiseless: *loc*, *causinf*; noisy: *change*, *wrm*; 4-6 distinct test functions per problem).

For the benchmark, we ran 50 independent runs of each algorithm on each test function, with randomized starting points and a budget of `500 * D` function evaluations. If an algorithm terminated before depleting the budget, it was restarted from a new random point.
We consider a run *successful* if the best achieved function value is within `eps` from the true optimum `f_opt` (or our best estimate thereof). We set a target of `eps` = 0.1 for noiseless functions and `eps` = 1 for noisy functions (smaller differences in log-likelihood are negligible due to data and model uncertainty).

## Results

**Fig 2: Performance on artificial test functions** ![Fig 2: Performance on artificial test functions](https://github.com/lacerbi/bads-dev/blob/master/docs/figures/fig2-demo.png "Fig 2: Performance on artificial test functions")

**A:** Noiseless *cec14* testbed. Plots of fraction of successful runs (`eps` < 0.1) as a function of number of function evaluations per number of dimensions, for different algorithms. Each panel contains a subset of test functions. *Left column*: Low-dimensional functions with `D = 2,3,5` (left). *Right column*: Medium-dimensional functions with `D = 10,20` (right). *Top row*: Smooth functions (sphere, ellipsoid, rotated ellipsoid and Rosenbrock). *Bottom row*: Non-smooth functions (step, Ackley, Griewank, Rastrigin). 
**B:** Fraction of successful runs (`eps` < 1) on the *cec14-noisy* testbed with additive Gaussian noise. BADS performs on par with or outperforms the best algorithms in each class (`fmincon` and MCS on smooth, noiseless problems; CMA-ES and MCS on non-smooth or noisy problems).

**Fig 3: Performance on behavioral model-fitting problems** ![Fig 3: Performance on behavioral model-fitting problems](https://github.com/lacerbi/bads-dev/blob/master/docs/figures/fig3-demo.png "Fig 3: Performance on behavioral model-fitting problems")

**A:** Maximization of noiseless log-likelihood functions of behavioral models; fraction of problems solved (see **Fig 2** for the legend). *Left panel*: Localization task modelled with location-dependent sensory noise (*loc*; 5 pairs of representative models and real datasets). *Right panel*: Perceptual causal inference task modelled with stimulus-dependent noise (*causalinf*; 5 models/real datasets). This problem is hard due to nested numerical integrations that yield a jagged log-likelihood landscape (for this we set a target tolerance of `eps` < 1, as per noisy problems).
**B:** Maximization of log-likelihood functions obtained through stochastic simulation (noisy); estimated fraction of problems solved after `500 * D` function evaluations. *Left panel*: Visual short-term memory change localization task (*change*; 5 synthetic datasets from meaningful parameter settings; `D = 2`). This problem is hard due to shallow log-likelihood landscapes with relatively high noise. *Right panel*: Bayesian model of a word-recognition memory task (*wrm*; 4 synthetic datasets; `D = 4`). 
**C:** Average execution time per function call. For most non-analytical models, the additional cost of running BADS is a small fraction of the model evaluation cost. The advantages of BADS are reduced when model evaluation is very quick (e.g., less than 0.1 s per function call, such as in the *loc* case), in which case it may be convenient to adopt a faster algorithm (although the choice of the alternative is not obvious).

## Conclusions

BADS [1] combines MADS [2] and Bayesian Optimization [3], similarly in principle to [4], but in a novel approach that makes it effective for model fitting. Previous work in Bayesian Optimization targets only slow simulations on much longer timescales [3,4]. 
BADS can be used out-of-the-box and will be soon available as a MATLAB package (a Python version may follow).

## References

1. Acerbi L & Ma WJ (2017). "Bayesian Adaptive Direct Search", in preparation.
2. Audet C & Dennis Jr JE (2006). "Mesh Adaptive Direct Search Algorithms for Constrained Optimization", *SIAM J Opt* 
([link](http://epubs.siam.org/doi/abs/10.1137/040603371)).
3. Snoek *et al*. (2012). "Practical Bayesian Optimization of Machine Learning Algorithms", *NIPS* 2012(4522) ([link](https://papers.nips.cc/paper/4522-practical-bayesian-optimization-of-machine-learning-algorithms)).
4. Gramacy RB & Le Digabel S (2015). "The mesh adaptive direct search algorithm with treed Gaussian process surrogates", *Pac J Opt* ([link](http://www.optimization-online.org/DB_HTML/2011/07/3090.html)).
5. Rios LM & Sahinidis NV (2013). "Derivative-free optimization: a review of algorithms and comparison of software implementations", *J Global Opt* ([link](http://link.springer.com/article/10.1007/s10898-012-9951-y)).

