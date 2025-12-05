# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

BADS (Bayesian Adaptive Direct Search) is a MATLAB optimization toolbox for solving difficult nonlinear optimization problems, particularly for model fitting with noisy or stochastic objective functions. It combines mesh adaptive direct search (MADS) with Gaussian process (GP) surrogate modeling.

## Commands

### Installation
```matlab
run install.m    % Add BADS to MATLAB path
```

### Testing
```matlab
bads('test')     % Run test suite (4 tests: deterministic, constrained, noisy, heteroskedastic)
```

### Running Examples
```matlab
bads_examples    % Interactive examples demonstrating various use cases
```

### Basic Usage
```matlab
[x, fval] = bads(@fun, x0, LB, UB, PLB, PUB);           % Basic optimization
[x, fval] = bads(@fun, x0, LB, UB, PLB, PUB, nonbcon);  % With non-bound constraints
[x, fval] = bads(@fun, x0, LB, UB, PLB, PUB, [], opts); % With options
```

### Releasing
```bash
./scripts/update_version.sh 1.2.0  # Updates version in README.md and bads.m
```
Then manually add a changelog entry at the end of `bads.m`.

## Architecture

BADS alternates between two main stages in its optimization loop:

1. **Poll Stage** (`poll/`): Evaluates points on a mesh around the current incumbent using directional basis vectors. Success expands the mesh; failure contracts it.
   - Default method: `pollMADS2N` - lower-triangular MADS with 2D random basis vectors

2. **Search Stage** (`search/`): Uses GP predictions to explore promising regions via acquisition function optimization.
   - Default method: `searchHedge` - probabilistically selects between evolution-strategy-inspired searches based on track record

### Key Components

- **`bads.m`**: Main entry point and optimization loop
- **`acq/`**: Acquisition functions for selecting evaluation points
  - Default: `acqLCB` (lower confidence bound) for both poll and search
- **`gpdef/`**: GP definition and hyperprior configuration (`gpdefBads`)
- **`gpml_fast/`**: Fast GP operations (training, prediction, updates)
- **`gpml-matlab-v3.6-2015-07-07/`**: Bundled GPML library for GP computations
- **`init/`**: Initialization functions (default: `initSobol`)
- **`utils/`**: Helper functions for bounds checking, coordinate transforms, etc.
- **`warp/`**: Warping functions (currently unsupported)

### Noise Handling

BADS supports three modes:
1. **Deterministic**: `options.UncertaintyHandling = false`
2. **Noisy (unknown variance)**: `options.UncertaintyHandling = true`, optionally set `options.NoiseSize`
3. **Heteroskedastic (known variance)**: `options.SpecifyTargetNoise = true` - objective returns `[fval, sd]`

## Key Options

```matlab
options = bads('defaults');              % Get default options struct
options.MaxFunEvals = 500*nvars;         % Function evaluation budget
options.MaxIter = 200*nvars;             % Iteration limit
options.TolMesh = 1e-6;                  % Mesh size termination tolerance
options.UncertaintyHandling = [];        % Auto-detect noise (true/false to specify)
options.PeriodicVars = [];               % Indices of periodic variables
options.Display = 'iter';                % 'iter', 'final', 'notify', or 'off'
```

## Testing Changes

After modifying core algorithms, run:
```matlab
bads('test')  % Should report 0 failed tests
```

For quick validation on a simple problem:
```matlab
[x, fval] = bads(@(x) sum(x.^2), zeros(1,3), -10*ones(1,3), 10*ones(1,3), -5*ones(1,3), 5*ones(1,3));
% Expected: x near [0,0,0], fval near 0
```
