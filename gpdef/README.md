# bads-dev
Bayesian Adaptive Direct Search - developer's version


### Installation instructions

- To install BADS, just clone or unpack the zipped repository where you want it and run the script `install.m`.
   - This will add the BADS base folder to the MATLAB search path.
   - No need to install other toolboxes.
   - BADS automatically handles the rest at runtime.
- To see if everything works, run `bads_test.m`.

### Documentation

- The BADS interface is similar to that of other MATLAB optimizers, such as `fmincon` and `patternsearch`.
- If you type `help bads` you will get usage information (to be extended).
   - You can also check the `bads_test.m` script.
- If you simply type `bads` you will get a default OPTIONS struct.
   - These are the OPTIONS you might *sometimes* want to play with (but in general you will be okay with the defaults). 
   - There are other hidden options which are not recommended for the user to modify.
- **BADS for noisy problems:** You need to set `OPTIONS.UncertaintyHandling = 1` and `OPTIONS.NoiseSize = sigma`. 
   - `sigma` is an estimate of the SD of the noise in your problem (it is not limited to this value).
   - The noise handling part still under testing and improvement.
