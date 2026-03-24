# Lie Groups

Lie group utilities for rotation representations used by attitude-aware
estimators. Currently provides SO(3) quaternion operations following the
Hamilton convention with w-first user-facing representation.

## Types

- [so3](so3.md) -- SO(3) rotation group: exp/log maps, compose, conjugate, normalize, skew

## When to use

Use the **so3** utilities when working with MEKF, manifold UKF, or any
algorithm that needs to map between rotation quaternions and their tangent
space (the Lie algebra so(3)). The exp and log maps provide the bridge
between 3-vector perturbations and unit quaternions.

## Theory

For the mathematical background, see
[Attitude Estimation Theory](../reference/attitude-estimation-theory.md).
