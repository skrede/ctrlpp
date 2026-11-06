# Lie Groups

Lie group utilities for rotation representations used by attitude-aware
estimators. Provides SO(3) quaternion operations following the Hamilton
convention with w-first user-facing representation; no other group is
implemented.

## Types

- [so3](so3.md)<br/> SO(3) rotation group: exp/log maps, compose, conjugate, normalize, skew, w-first serialization

## When to use

Use the **so3** utilities when working with MEKF, manifold UKF, or any
algorithm that needs to map between rotation quaternions and their tangent
space (the Lie algebra so(3)). The exp and log maps provide the bridge
between 3-vector perturbations and unit quaternions.

## Theory

For the mathematical background, see
[SO(3) and Lie Groups](../../background/so3.md) and
[Attitude Estimation](../../background/attitude-estimation.md).
