# SPHKernels.jl

This package supplies a number of kernels frequently used in Smoothed-Particle Hydrodynamics (SPH), as well as functions to evaluate their values and derivatives in 2D and 3D.

The implementation closely follows the one in [Gadget2](https://wwwmpa.mpa-garching.mpg.de/gadget/), see [Springel (2005)](https://ui.adsabs.harvard.edu/abs/2005MNRAS.364.1105S/abstract) for details.

These kernels include the B-splines (`Cubic` and `Quintic`) suggested in [Monaghan & Lattanzio (1985)](https://ui.adsabs.harvard.edu/abs/1985A%26A...149..135M/abstract), the Wendland functions (`WendlandC2`, `WendlandC4` and `WendlandC6` from [Wendland (2009)](https://www.researchgate.net/publication/220179293_Divergence-Free_Kernel_Methods_for_Approximating_the_Stokes_Problem)) as suggested in [Dehnen & Aly (2012)](https://academic.oup.com/mnras/article/425/2/1068/1187211) and `WendlandC8` as suggested by [Kummer et. al. (2019)](https://arxiv.org/abs/1902.02330).


> :warning: **The version numbering of this package is unfortunately not really reflective of the state. I made an error on the original setup of the repository, so I had to start out with version 1.0. View this more as v0.2, instead of v2.0!**: Please sanity-check everything before you use it in production!
     

# Table of contents

```@meta
CurrentModule = SPHKernels
DocTestSetup = quote
    using SPHKernels
end
```

```@contents
Pages = [ "index.md", 
          "install.md", 
          "kernels.md",
          "extending.md", 
          "api.md"]
Depth = 3
```