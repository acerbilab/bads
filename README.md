# Bayesian Adaptive Direct Search (BADS)

## What is it

BADS is a novel, fast Bayesian optimization algorithm for MATLAB designed to solve difficult optimization problems, in particular related to fitting computational models (e.g., via [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation)).

In our benchmark on real model-fitting problems, BADS performed on par or better than 15 other common and state-of-the-art MATLAB optimizers, such as `fminsearch`, `fmincon`, and `cmaes` [1].

BADS is recommended when no gradient information is available, and the objective function is non-analytical or *noisy*, for example evaluated through numerical approximation or via simulation. 

BADS requires no specific tuning and runs off-the-shelf like other built-in MATLAB optimizers such as `fminsearch`.

## Installation

[**Download the latest version of BADS as ZIP file**](https://github.com/lacerbi/bads/archive/master.zip).
- To install BADS, clone or unpack the zipped repository where you want it and run the script `install.m`.
   - This will add the BADS base folder to the MATLAB search path.
- To see if everything works, run `bads('test')`.

## Quick start

The BADS interface is similar to that of other MATLAB optimizers. The basic usage is:

```matlab
[X,FVAL] = bads(FUN,X0,LB,UB,PLB,PUB);
```
with input parameters:
- `FUN`, a function handle to the objective function to minimize (typically, the log likelihood of a dataset and model, for a given input parameter vector);
- `X0`, the starting point of the optimization (a row vector);
- `LB` and `UB`, hard lower and upper bounds (can contain `-Inf` and `Inf`);
- `PLB` and `PUB`, *plausible* lower and upper bounds, that is a box where you would expect to find almost all solutions.

The output parameters are:
- `X`, the found optimum.
- `FVAL`, the (estimated) function value at the optimum.

For more usage examples, see [**bads_examples.m**](https://github.com/lacerbi/bads-dev/blob/master/bads_examples.m). You can also type `help bads` to see the documentation.

## How does it work

BO approximates the target function via a Gaussian Process model, and searches the parameter landscape with a trade-off between exploration of uncertain regions and exploitation of promising solutions.

This project is under active development. If you find a bug, or anything that needs correction, please let me know.


### Reference

1. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian Optimization for Model Fitting with Bayesian Adaptive Direct Search. *arXiv preprint*, arXiv:YYMM.NNNN. (link)

You can cite BADS in your work with something along the lines of

> We optimized the log likelihoods of our models using Bayesian adaptive direct search (BADS; Acerbi and Ma, 2017). BADS Briefly, BADS alternates between a series of fast, local BO steps (the \search{} stage of MADS) and, when the search repeatedly fails, a systematic, slower exploration of the mesh grid (\poll{} stage). 

If you use BADS, you can demonstrate your appreciation by

- *Starring* the BADS repository on GitHub;
- Follow me on Twitter for more updates about BADS and other projects I am involved;
- Cite BADS as

