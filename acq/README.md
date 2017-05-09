# Acquisition functions

The *acquisition function* provides an estimate of the (negative) value of evaluating a point for the purpose of optimization.
BADS selects the next point to evaluate, both during POLL and during SEARCH, by minimizing the chosen acquisition function (POLL and SEARCH can use different acquisition functions). 

Typically, the evaluation function takes into account both the mean and the uncertainty associated with the current  Gaussian process approximation of the objective function, setting a trade-off between exploitation and exploration.

The default acquisition function used in BADS for both POLL and SEARCH is `acqLCB`.
Alternative acquisition functions can be set in BADS via `OPTIONS.PollAcqFcn` and `OPTIONS.SearchAcqFcn`. 
For the standard user, there is no reason to try out other acquisition functions.

## List of supported acquisition functions

- `acqLCB`: lower confidence bound (LCB), default [1]. Takes as additional argument the value of `sqrtbeta`, which otherwise is computed according to Eq. 2 in the paper. For example usage:
```matlab
sqrtbeta = [];                               % Could be a scalar
OPTIONS.SearchAcqFcn = {@acqLCB,sqrtbeta};   % Set acquisition function for SEARCH
```
- `acqNegEI`: negative expected improvement (EI) [2]. The common choice for Bayesian optimization, but not the best-performing metric in BADS.
- `acqNegPI`: negative probability of improvement (PI) [3]. Generally not recommended, too greedy.
- `acqRnd`: random noise (provides no information). Used for debugging purposes.

Other acquisition functions are either under development or currently unsupported.

### References

1. Srinivas, N., Krause, A., Kakade, S. M., & Seeger, M. (2009). Gaussian process optimization in the bandit setting: No regret and experimental design. *Proceedings of the 27th International Conference on Machine Learning (ICML-10)*, pp. 1015-1022. ([link](http://machinelearning.wustl.edu/mlpapers/paper_files/icml2010_SrinivasKKS10.pdf))
2. Mockus, J.,  Tiesis, V., & Zilinskas, A. (1978). *Toward Global Optimization*, volume 2, chapter The Application of Bayesian Methods for Seeking the Extremum, pp. 117–128. Elsevier.
3. Kushner, H. J. (1964). A new method of locating the maximum of an arbitrary multipeak curve in the presence of noise. *J. Basic Engineering*, **86**:97–106. ([link](http://fluidsengineering.asmedigitalcollection.asme.org/article.aspx?articleid=1431594))
