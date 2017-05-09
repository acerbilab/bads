# Warping functions

An (output) *warping function* is a nonlinear transformation applied to the observed function values, generally to restrict the range of observed outputs and to obtain a better GP approximation (e.g. to model non-Gassian noise) [1]. In Evolution Strategies, a similar concept is called *fitness shaping* [2]. Note that this differs from *input warping*, another fairly common strategy in Bayesian optimization [3].

BADS does not currently support warping functions, as an early exploration of this feature did not show a consistent improvement over vanilla GPs. Unless there are new developments, this folder may eventually be removed in future versions of BADS.

### References:

1) Snelson, E., Rasmussen, C. E., & Ghahramani, Z. (2003). Warped Gaussian Processes. In *NIPS* (pp. 337-344). ([link](https://papers.nips.cc/paper/2481-warped-gaussian-processes.pdf))
2) Wierstra, D., Schaul, T., Peters, J., & Schmidhuber, J. (2008). Natural evolution strategies. In *Evolutionary Computation, 2008. CEC 2008.(IEEE World Congress on Computational Intelligence). IEEE Congress on* (pp. 3381-3387). IEEE. ([link](http://ieeexplore.ieee.org/abstract/document/4631255/))
3) Snoek, J., Swersky, K., Zemel, R. S., & Adams, R. P. (2014). Input Warping for Bayesian Optimization of Non-Stationary Functions. In *ICML* (pp. 1674-1682). ([link](http://proceedings.mlr.press/v32/snoek14.pdf))
