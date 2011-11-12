# Gibbs Sampler in R using proto

This contains R code to implement Gibbs samplers using a
prototype-based programming paradigm with the **proto** package.  See
the vignette in the **proto** package for a general discussion of
prototype based programming and how to use proto objects.

`GibbsSampler` in `protoGibbs.R` provides an object from which other
samplers can be cloned. See `normalLinear.R` and `gibbsDLM.R` for two
examples of Gibbs samplers.

This implementation was inspired by the now unmaintained Umacs R
package.

