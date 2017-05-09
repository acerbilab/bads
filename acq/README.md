# Acquisition functions

The *acquisition function* provides an estimate of the (negative) value of evaluating a point for the purpose of optimization.
BADS selects the next point to evaluate, both during POLL and during SEARCH, by minimizing the chosen acquisition function (POLL and SEARCH can use different acquisition functions). 

Typically, the evaluation function takes into account both the mean and the uncertainty associated with the current  Gaussian process approximation of the objective function, to counterbalance exploitation and exploration.

The default acquisition function used in BADS for both POLL and SEARCH is `acqLCB`.
Alternative acquisition functions can be set in BADS via `OPTIONS.PollAcqFcn` and `OPTIONS.SearchAcqFcn`. 
For the standard user, there is no reason to try out other acquisition functions.

## List of supported acquisition functions

- `acqLCB`: lower confidence bound (LCB), default. Takes as additional argument the value of `sqrtbeta`, which otherwise is computed according to Eq. 3 in the paper. For example usage:
```matlab
sqrtbeta = [];                               % Could be a scalar
OPTIONS.SearchAcqFcn = {@acqLCB,sqrtbeta};   % Set acquisition function for SEARCH
```
- `acqNegEI`: negative expected improvement (EI). The common choice for Bayesian optimization, but not the best-performing metric in BADS.
- `acqNegPI`: negative probability of improvement (PI). Generally not recommended, too greedy.
- `acqRnd`: random noise (provides no information). Used for debugging purposes.

Other acquisition functions are either under development or currently unsupported.
