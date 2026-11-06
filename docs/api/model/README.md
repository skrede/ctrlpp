# Model Utilities

Types and concepts for representing, converting, analyzing, and propagating
dynamic system models. These are the building blocks that controllers,
estimators, and identification algorithms operate on.

## Types

### Representations

- [state_space](state-space.md)<br/> Linear continuous-time and discrete-time state-space models (A, B, C, D matrices)
- [transfer_function](transfer-function.md)<br/> SISO transfer function representation (numerator/denominator polynomials)

### Operations

- [discretize](discretize.md)<br/> Continuous-to-discrete conversion (ZOH, Tustin with optional prewarping, forward and backward Euler)
- [conversion](conversion.md)<br/> SISO transfer function to state-space and back, with `conversion_error` refusals
- [analysis](analysis.md)<br/> Stability, controllability, and observability checks
- [propagate](propagate.md)<br/> State propagation and output computation utilities

### Concepts

These are C++20 concepts. `dynamics_model` and `measurement_model` constrain the
model parameters of the EKF, UKF, particle filter, MHE and NMHE (NMPC requires
`dynamics_model` only); the differentiable refinements select analytic Jacobians
over finite differences. They define what a user-provided model must look like.

- [dynamics_model](dynamics-model.md)<br/> Callable `f(x, u)` returning next state
- [measurement_model](measurement-model.md)<br/> Callable `h(x)` returning predicted measurement
- [differentiable_dynamics](differentiable-dynamics.md)<br/> Dynamics with analytical Jacobian
- [differentiable_measurement](differentiable-measurement.md)<br/> Measurement with analytical Jacobian
- [constraint_model](constraint-model.md)<br/> Shape of the path and terminal constraint callables NMPC takes

## When to use

Use **state_space** and **transfer_function** for linear system representations.
Convert between them with **conversion**, discretize continuous models with
**discretize**, and check properties with **analysis**.

Use **propagate** when you need explicit state propagation outside of a filter
or controller (e.g., in simulation).

Implement the **concepts** when writing custom dynamics for EKF, UKF, NMPC, or
MHE &mdash; the compiler checks your model satisfies the required interface at
compile time.

The convenience header `<ctrlpp/model.h>` includes every header of the module.
`<ctrlpp/analysis.h>`, `<ctrlpp/conversion.h>`, `<ctrlpp/discretize.h>`,
`<ctrlpp/propagate.h>`, `<ctrlpp/state_space.h>` and `<ctrlpp/transfer_function.h>`
each forward one of them.
