# Bayesian Adaptive Direct Search (BADS) - v1.0.3

## What is it

BADS is a novel, fast Bayesian optimization algorithm designed to solve difficult optimization problems, in particular related to fitting computational models (e.g., via [maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation)).

BADS has been intensively tested for fitting behavioral, cognitive, and neural models and is currently being used in more than a dozen projects in the [Ma lab](http://www.cns.nyu.edu/malab/).
In our benchmark with real model-fitting problems, BADS performed on par or better than many other common and state-of-the-art MATLAB optimizers, such as `fminsearch`, `fmincon`, and `cmaes` [[1](#reference)].

BADS is recommended when no gradient information is available, and the objective function is non-analytical or *noisy*, for example evaluated through numerical approximation or via simulation. 

BADS requires no specific tuning and runs off-the-shelf like other built-in MATLAB optimizers such as `fminsearch`.


## Installation

[**Download the latest version of BADS as a ZIP file**](https://github.com/lacerbi/bads/archive/master.zip).
- To install BADS, clone or unpack the zipped repository where you want it and run the script `install.m`.
   - This will add the BADS base folder to the MATLAB search path.
- To see if everything works, run `bads('test')`.

## Quick start

The BADS interface is similar to that of other MATLAB optimizers. The basic usage is:

```matlab
[X,FVAL] = bads(FUN,X0,LB,UB,PLB,PUB);
```
with input parameters:
- `FUN`, a function handle to the objective function to minimize (typically, the *negative* log likelihood of a dataset and model, for a given input parameter vector);
- `X0`, the starting point of the optimization (a row vector);
- `LB` and `UB`, hard lower and upper bounds;
- `PLB` and `PUB`, *plausible* lower and upper bounds, that is a box where you would expect to find almost all solutions.

The output parameters are:
- `X`, the found optimum.
- `FVAL`, the (estimated) function value at the optimum.

For more usage examples, see [**bads_examples.m**](https://github.com/lacerbi/bads/blob/master/bads_examples.m). You can also type `help bads` to display the documentation.

For practical recommendations, such as how to set `LB` and `UB`, and any other question, check out the FAQ on the [BADS wiki](https://github.com/lacerbi/bads/wiki).

## How does it work

BADS follows a [mesh adaptive direct search](http://epubs.siam.org/doi/abs/10.1137/040603371) (MADS) procedure for function minimization that alternates **poll** steps and **search** steps (see **Fig 1**). 

- In the **poll** stage, points are evaluated on a mesh by taking steps in one direction at a time, until an improvement is found or all directions have been tried. The step size is doubled in case of success, halved otherwise. 
- In the **search** stage, a [Gaussian process](https://en.wikipedia.org/wiki/Gaussian_process) (GP) is fit to a (local) subset of the points evaluated so far. Then, we iteratively choose points to evaluate according to a *lower confidence bound* strategy that trades off between exploration of uncertain regions (high GP uncertainty) and exploitation of promising solutions (low GP mean).

**Fig 1: BADS procedure** ![BADS procedure](https://github.com/lacerbi/bads/blob/master/docs/bads-cartoon.png "Fig 1: BADS procedure")

See [here](https://github.com/lacerbi/optimviz) for a visualization of several optimizers at work, including BADS.

See our paper for more details [[1](#reference)].

## Troubleshooting

If you have trouble doing something with BADS:

- Check out the FAQ on the [BADS wiki](https://github.com/lacerbi/bads/wiki);
- Contact me at <luigi.acerbi@nyu.edu>, putting 'BADS' in the subject of the email.

This project is under active development. If you find a bug, or anything that needs correction, please let me know.

## Reference

1. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian Optimization for Model Fitting with Bayesian Adaptive Direct Search. *arXiv preprint*, arXiv:1705.04405. ([link](https://arxiv.org/abs/1705.04405))

You can cite BADS in your work with something along the lines of

> We optimized the log likelihoods of our models using Bayesian adaptive direct search (BADS; Acerbi and Ma, 2017). BADS alternates between a series of fast, local Bayesian optimization steps and a systematic, slower exploration of a mesh grid. 

Besides formal citations, you can demonstrate your appreciation for BADS in the following ways:

- *Star* the BADS repository on GitHub;
- [Follow me on Twitter](https://twitter.com/AcerbiLuigi) for updates about BADS and other projects I am involved with;
- Tell me about your model-fitting problem and your experience with BADS (positive or negative) at <luigi.acerbi@nyu.edu> (putting  'BADS' in the subject of the email).

### License

BADS is released under the terms of the [GNU General Public License v3.0](https://github.com/lacerbi/bads/blob/master/LICENSE.txt).
