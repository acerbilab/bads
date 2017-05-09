# Poll generating functions

The *poll generating function* determines how the polling set is built during the POLL stage.

The default poll generating function is set as follows

```
OPTIONS.PollMethod = @pollMADS2N;
```

## List of poll generating functions

- `pollMADS2N`: Lower-triangular MADS poll generation function (default), with `2*D` random basis vectors.
- `pollGPS2N`: Generalized pattern search axis-aligned poll generation function, with `2*D` basis vectors.

### References

1. Audet, C., & Dennis Jr, J. E. (2006). Mesh adaptive direct search algorithms for constrained optimization. *SIAM Journal on Optimization*, **17**(1), 188-217. ([link](http://www.caam.rice.edu/caam/trs/2004/TR04-02.pdf))
2. Kolda, T. G., Lewis, R. M., & Torczon, V. (2003). Optimization by direct search: New perspectives on some classical and modern methods. *SIAM Review*, **45**(3), 385-482. ([link](http://www.cs.wm.edu/~va/research/sirev.pdf))
