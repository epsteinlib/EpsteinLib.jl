# `EpsteinZetaFunction.jl`

This is a julia interface for [https://github.com/epsteinlib/epsteinlib](https://github.com/epsteinlib/epsteinlib).

Currently waiting for the wrapper `Epsteinlib_jll.jl`. This will work in all architectures except `aarch64-apple-darwin` due to a linking error in `BinaryBuilder.jl`.

For `aarch64-apple-darwin` the library is hand-compiled from `epsteinlib` using commit `2a367d4`.
