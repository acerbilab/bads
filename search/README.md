# Search functions

The *search function* determines how the search is conducted in the SEARCH stage.
The default in BADS is a *hedging* search, `searchHedge`, that chooses between two kinds of *evolution strategies* inspired search, `searchES` (type 1 and type 2).

```matlab
OPTIONS.SearchMethod = {@searchHedge,{{@searchES,1,1},{@searchES,2,1}}};
```

The default SEARCH method has been thoroughly tested, so unless there is a specific rationale we do not recommend to use alternative search functions.

## List of search functions

- `searchHedge`: Meta-search function that chooses probabilistically between other search functions, based on their track record of cumulative improvement, according to the **Hedge** algorithm [1] (default).
- `searchES`: [Evolution strategies](http://www.scholarpedia.org/article/Evolution_strategies) inspired search, with different subtypes that differ in the way the search covariance matrix is constructed; such as (1) weighted covariance matrix [2]; (2) GP-rescaled diagonal matrix; (3) identity matrix.
- `searchGauss`: A simple isotropic multivariate normal random search around the current incumbent.
- `searchWCM`: A single weighted covariance matrix adaptation search step (inspired by CMA-ES [2]).
- `searchCombine`: Meta-search function that combines multiple search functions, dividing samples among them.
- `searchOptim`: Local optimization of expected improvement via `fmincon` (deprecated).
- `searchGrid`: Build space-filling quasi-random grid surrounding the current incumbent.

Other unlisted search functions are not currently supported.

### References

1. Hoffman, M. D., Brochu, E., & de Freitas, N. (2011). Portfolio Allocation for Bayesian Optimization. In *UAI* (pp. 327-336). ([link](https://pdfs.semanticscholar.org/1a7f/d7b566697c9b69e64b27b68db4384314d925.pdf))
2. Hansen, N., Müller, S. D., & Koumoutsakos, P. (2003). Reducing the time complexity of the derandomized evolution strategy with covariance matrix adaptation (CMA-ES). *Evolutionary Computation*, **11**(1), 1-18. ([link](https://www.lri.fr/~hansen/evco_11_1_1_0.pdf))
