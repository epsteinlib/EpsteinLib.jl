# SPDX-FileCopyrightText: 2025 Jonathan Busse <jonathan@jbusse.de>
# SPDX-License-Identifier: AGPL-3.0-only
#
# @file test_epsteinzeta.jl
# @brief Simple numerical test for the Epstein zeta function wrapper.
#
# This script loads `EpsteinZetaFunction.jl`, evaluates the Epstein zeta function
# for a fixed set of parameters, compares the result to a known reference value,
# and exits with code 0 on success (|Δ| < 1e-14) or 1 on failure.
#
# Run with:
#     julia test_epsteinzeta.jl

include("EpsteinZetaFunction.jl")
using .EpsteinZetaFunction
using LinearAlgebra
using Printf

ν = 1.5
A = [1.0 0.1;
     0.0 2.0]
x = [0.0, 0.0]
y = [0.5, 0.5]

ref = -1.5212136988975367
z = epsteinzeta(ν, A, x, y)
val = real(z)

@printf("Epstein zeta value:       %.16e\n", val)
@printf("Reference value:          %.16e\n", ref)
@printf("Absolute difference:      %.2e\n", abs(val - ref))

tol = 1e-14
success = abs(val - ref) < tol

println(success ? "Test passed ✅" : "Test failed ❌ (|Δ| = $(abs(val - ref)))")
exit(success ? 0 : 1)
