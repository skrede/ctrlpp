# Model Utilities

Types and concepts for representing, converting, analysing, and propagating
dynamic system models. These are the building blocks that controllers,
estimators, and identification algorithms operate on.

## Types

### Representations

- [state_space](state-space.md) -- Linear discrete state-space model (A, B, C, D matrices)
- [transfer_function](transfer-function.md) -- Transfer function representation (numerator/denominator polynomials)

### Operations

- [discretise](discretise.md) -- Continuous-to-discrete conversion (ZOH via Van Loan matrix exponential)
- [conversion](conversion.md) -- Transfer function to state-space and back
- [analysis](analysis.md) -- Stability, controllability, and observability checks
- [propagate](propagate.md) -- State propagation and output computation utilities

### Concepts

These are C++23 concepts used as template constraints by EKF, UKF, particle
filter, NMPC, MHE, and other algorithm types. They define what a user-provided
model must look like.

- [dynamics_model](dynamics-model.md) -- Callable `f(x, u)` returning next state
- [measurement_model](measurement-model.md) -- Callable `h(x)` returning predicted measurement
- [differentiable_dynamics](differentiable-dynamics.md) -- Dynamics with analytical Jacobian
- [differentiable_measurement](differentiable-measurement.md) -- Measurement with analytical Jacobian
- [constraint_model](constraint-model.md) -- Path and terminal constraint callables for NMPC

## When to use

Use **state_space** and **transfer_function** for linear system representations.
Convert between them with **conversion**, discretise continuous models with
**discretise**, and check properties with **analysis**.

Use **propagate** when you need explicit state propagation outside of a filter
or controller (e.g., in simulation).

Implement the **concepts** when writing custom dynamics for EKF, UKF, NMPC, or
MHE -- the compiler checks your model satisfies the required interface at
compile time.
