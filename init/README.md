# Initial design

The *initial design* function specifies how points are selected in the initialization stage. These points are used, together with the provided `x0`, as a training set for the Gaussian process model.

## List of initial design functions

- `initSobol`: initialization via a space-filling quasi-random [Sobol sequence](https://en.wikipedia.org/wiki/Sobol_sequence) in the plausible box specified by `PLB` and `PUB`. This is the default option.
- `initLHS`: initialization via a [Latin hypercube sampling](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) in the plausible box.
- `initRand`: initialization via a random uniform sampling in the plausible box.
