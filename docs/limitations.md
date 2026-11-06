# Limitations

A limitation entry states something an API does not do, or does only under a condition, that a reader could
otherwise take for granted. Each entry lives on the page of the API it limits, under that page's Limitations
heading or beside the claim it qualifies, and this page links every one of them, grouped by module.

## Control

- [pid: `set_params` needs a copy-assignable callable](api/control/pid/feed-forward.md#set_params-needs-a-copy-assignable-callable) -- `set_params` copy-assigns the configuration, so it does not compile with a capturing-lambda feed-forward.
- [pid: a non-finite metric refuses the control cycle](api/control/pid/performance.md#a-non-finite-metric-refuses-the-control-cycle) -- a performance metric that overflows refuses the cycle with `non_finite_result` although the command was finite.
- [mrac: the adaptation fixes the Lyapunov weight to the identity](api/control/mrac.md#the-adaptation-fixes-the-lyapunov-weight-to-the-identity) -- the update uses `P = I`, a valid weight for the discrete reference model only when `I - A_m^T A_m` is positive definite; not configurable, not checked.
- [mrac: the update carries no sample-time factor](api/control/mrac.md#the-update-carries-no-sample-time-factor) -- the adaptation gains absorb the sample time, and the update is not normalized by signal magnitude.
- [mrac: MIMO adaptation assumes the reference input matrix is the plant's up to a positive scalar](api/control/mrac.md#mimo-adaptation-assumes-the-reference-input-matrix-is-the-plants-up-to-a-positive-scalar) -- with several inputs the decrease argument holds only under that assumption.
- [l1: the adaptation law fixes the Lyapunov weight to the identity](api/control/l1.md#the-adaptation-law-fixes-the-lyapunov-weight-to-the-identity) -- the same `P = I` specialization as `mrac`, with no sample-time factor.
- [l1: the adaptation gain has a ceiling](api/control/l1.md#the-adaptation-gain-has-a-ceiling) -- the discrete update is stable only for `0 < gamma_i b_i^2 < 2 (1 + a_i)`, unlike the continuous-time law.
- [l1: the filter condition of the L1 theory is not checked](api/control/l1.md#the-filter-condition-of-the-l1-theory-is-not-checked) -- the L1 performance bounds need a filter and uncertainty condition the controller cannot test.
- [place: the gain is not checked against the requested poles](api/control/place.md#the-gain-is-not-checked-against-the-requested-poles) -- an ill-conditioned controllability matrix passes the rank test, and the returned gain can miss the poles.
- [dare: a singular state matrix is refused](api/control/dare.md#a-singular-state-matrix-is-refused) -- a plant with an eigenvalue at zero is refused, and so are the LQR, LQI and MPC terminal-cost designs on it.
- [lqr: the finite-horizon recursions do not test their factorization](api/control/lqr.md#the-finite-horizon-recursions-do-not-test-their-factorization) -- `lqr_finite` and `lqr_tv_gains` solve a possibly singular `R + B'PB` without a rank test.
- [lqr: the sign-function tag's `warmup_iterations` is not forwarded](api/control/lqr.md#the-sign-function-tags-warmup_iterations-is-not-forwarded) -- `care` and `lqr_gain_continuous` run the sign-function solve with `warmup_iterations = 3` whatever the tag says.

## Estimation

- [kalman: the feedthrough matrix is not used](api/estimation/kalman.md#the-feedthrough-matrix-is-not-used) -- the predicted measurement is `C x`; subtract `D u` before `update`.
- [luenberger: the feedthrough matrix is not used](api/estimation/luenberger.md#the-feedthrough-matrix-is-not-used) -- the correction compares the measurement with `C x`.
- [ukf: the Merwe defaults are not tuned for dimension, scale or precision](api/estimation/ukf.md#the-merwe-defaults-are-not-tuned-for-dimension-scale-or-precision) -- the default sigma-point parameters are fixed constants, which degrade `float` filters and large state magnitudes.
- [particle_filter: the default likelihood assumes additive Gaussian measurement noise](api/estimation/particle-filter.md#the-default-likelihood-assumes-additive-gaussian-measurement-noise) -- non-Gaussian or wrapped measurements need a likelihood policy; the process noise is always additive Gaussian.
- [particle_filter: `map_estimate()` returns the first particle after a resample](api/estimation/particle-filter.md#map_estimate-returns-the-first-particle-after-a-resample) -- with every weight equal, the estimate is chosen by array position.
- [Extended Kalman filter theory: limitations of the method](background/ekf.md#limitations) -- first-order accuracy, the Jacobian requirement, no general convergence guarantee and the Gaussian assumption.

## Predictive control and moving-horizon estimation

- [mpc: no offset-free formulation exists](api/mpc/mpc.md#no-offset-free-formulation-exists) -- there is no disturbance model or estimator coupling, so a constant disturbance leaves a steady-state offset.
- [mpc: tracking has no steady-state input target](api/mpc/mpc.md#tracking-has-no-steady-state-input-target) -- the input is penalized about zero, so a reference that needs a nonzero input settles short of it.
- [mpc: the softened state bounds are not an exact penalty](api/mpc/mpc.md#the-softened-state-bounds-are-not-an-exact-penalty) -- the quadratic slack penalty violates an active bound even when the hard problem is feasible.
- [mpc: each terminal-set function computes a different set](api/mpc/mpc.md#each-terminal-set-function-computes-a-different-set) -- one routine computes a robust control-invariant set and the other the maximal output-admissible set, both for at most four states.
- [mpc: terminal sets](api/mpc/mpc.md#terminal-sets) -- the inscribed box is not invariant, and the ellipsoid's state faces are not checked to straddle zero.
- [MPC theory: the terminal constraint set](background/mpc.md#terminal-constraint-set-mathcalx_f) -- what each terminal-set construction guarantees.
- [nmpc: no offset-free formulation exists](api/mpc/nmpc.md#no-offset-free-formulation-exists) -- as for `mpc`, a constant disturbance leaves a steady-state offset.
- [nmpc: tracking has no steady-state input target](api/mpc/nmpc.md#tracking-has-no-steady-state-input-target) -- the inputs are weighted about zero.
- [nmpc: state and input bounds are never softened](api/mpc/nmpc.md#state-and-input-bounds-are-never-softened) -- softening reaches only the path and terminal constraints.
- [nmpc: the two backends' augmented-Lagrangian options are not interchangeable](api/mpc/nmpc.md#the-two-backends-augmented-lagrangian-options-are-not-interchangeable) -- NLopt's variant penalizes the equalities only; argmin's folds the inequalities in as well.
- [nmpc: the compile-time horizon](api/mpc/nmpc.md#nmpc-compile-time-horizon) -- the compile-time-horizon controller adds no slack variables, so a path or terminal constraint needs `soft_constraints = false`.
- [mhe: the window shares one linearization point](api/mpc/mhe.md#the-window-shares-one-linearization-point) -- the model is linearized once per solve and every stage shares it; `nmhe` solves the nonlinear problem.
- [mhe: the linearized model drops the affine offset](api/mpc/mhe.md#the-linearized-model-drops-the-affine-offset) -- a model with a constant offset is estimated with the wrong model.
- [mhe: the solver cannot be configured](api/mpc/mhe.md#the-solver-cannot-be-configured) -- the solver is default-constructed; a setting needs a wrapper type.
- [nmhe: `create` does not report a solver setup failure](api/mpc/nmhe.md#create-does-not-report-a-solver-setup-failure) -- a setup failure surfaces only as a failed status on every window solve.
- [nmhe: the solver cannot be configured](api/mpc/nmhe.md#the-solver-cannot-be-configured) -- as for `mhe`.
- [nmhe: only the path constraint is softened](api/mpc/nmhe.md#only-the-path-constraint-is-softened) -- the state bounds and the residual bound stay hard.
- [nmhe: the NLP derivatives are always finite differences](api/mpc/nmhe.md#the-nlp-derivatives-are-always-finite-differences) -- the model Jacobians are used only by the embedded EKF.
- [osqp_solver: the adapter checks the bound lengths, not the cost vector or the warm start](api/mpc/osqp-solver.md#the-adapter-checks-the-bound-lengths-not-the-cost-vector-or-the-warm-start) -- the other vectors are read for the posed length unchecked.
- [argmin_qp_solver: the adapter leaves every dimension check to argmin](api/mpc/argmin-qp-solver.md#the-adapter-leaves-every-dimension-check-to-argmin) -- and a refused re-setup leaves the earlier problem posed.
- [nlopt_solver: the auglag variants penalize only the equality rows](api/mpc/nlopt-solver.md#the-auglag-variants-penalize-only-the-equality-rows) -- the inequalities and bounds go to the inner solver.
- [nlopt_solver: the adapter checks no dimension itself](api/mpc/nlopt-solver.md#the-adapter-checks-no-dimension-itself) -- the constraint row bounds are read unchecked.
- [nlopt_solver: the constraints are differentiated by central differences](api/mpc/nlopt-solver.md#the-constraints-are-differentiated-by-central-differences) -- a supplied constraint Jacobian is not read.
- [nlopt_solver: only SLSQP solved the shipped NMPC transcription at the default settings](api/mpc/nlopt-solver.md#only-slsqp-solved-the-shipped-nmpc-transcription-at-the-default-settings) -- the measured outcome of every other algorithm.
- [argmin_solver: `argmin_auglag` penalizes the inequality rows as well as the equalities](api/mpc/argmin-solver.md#argmin_auglag-penalizes-the-inequality-rows-as-well-as-the-equalities) -- not the same method as NLopt's.
- [argmin_solver: a bound-only policy drops the constraint rows](api/mpc/argmin-solver.md#a-bound-only-policy-drops-the-constraint-rows) -- the solve ignores the rows and reports no violation.
- [argmin_solver: the adapter checks dimensions only at setup](api/mpc/argmin-solver.md#the-adapter-checks-dimensions-only-at-setup) -- a solve before the first setup is undefined behavior.

## System identification

- [batch_arx: the record checks are not an identifiability test](api/sysid/batch-arx.md#the-record-checks-are-not-an-identifiability-test) -- the refusals only make the arithmetic defined; rank is reported, never refused.

## Model utilities

- [analysis: a false answer means not shown to hold](api/model/analysis.md#a-false-answer-means-not-shown-to-hold) -- a non-finite input or an eigensolver that did not converge also reads as `false`.
- [analysis: rank is decided within rounding](api/model/analysis.md#rank-is-decided-within-rounding) -- a system rank-deficient within the pivot threshold reads as uncontrollable or unobservable.
- [analysis: the stability tests have no margin](api/model/analysis.md#the-stability-tests-have-no-margin) -- the pole comparisons are exact.

## Trajectory

- [trapezoidal_trajectory: a boundary velocity above `v_max` is refused, with zero tolerance by default](api/trajectory/trapezoidal-trajectory.md#a-boundary-velocity-above-v_max-is-refused-with-zero-tolerance-by-default) -- the allowance `boundary_velocity_rounding_ops` defaults to zero.
- [trapezoidal_trajectory: composing segments through via points is not supported](api/trajectory/trapezoidal-trajectory.md#composing-segments-through-via-points-is-not-supported) -- the library computes no join velocity.
- [trapezoidal_trajectory: a zero displacement with unequal boundary speeds summing to a positive value is refused](api/trajectory/trapezoidal-trajectory.md#a-zero-displacement-with-unequal-boundary-speeds-summing-to-a-positive-value-is-refused) -- under either policy.
- [trapezoidal_trajectory: the velocity limit is a precondition](api/trajectory/trapezoidal-trajectory.md#the-velocity-limit-is-a-precondition)
- [double_s_trajectory: a boundary velocity above `v_max` is refused, with zero tolerance by default](api/trajectory/double-s-trajectory.md#a-boundary-velocity-above-v_max-is-refused-with-zero-tolerance-by-default) -- as for the trapezoidal profile.
- [double_s_trajectory: composing segments through via points is not supported](api/trajectory/double-s-trajectory.md#composing-segments-through-via-points-is-not-supported) -- as for the trapezoidal profile.
- [double_s_trajectory: a zero displacement is realized only from rest to rest](api/trajectory/double-s-trajectory.md#a-zero-displacement-is-realized-only-from-rest-to-rest)
- [double_s_trajectory: the velocity limit is a precondition](api/trajectory/double-s-trajectory.md#the-velocity-limit-is-a-precondition)
- [Multi-waypoint paths: chaining point-to-point profiles through via points is not supported](guides/trajectory/multi-waypoint.md#chaining-point-to-point-profiles-through-via-points-is-not-supported) -- the splines, the B-spline and `piecewise_trajectory` are the supported routes.

## Error reporting

- [`ctrlpp::expected` implements part of the `std::expected` interface](guides/patterns/error-reporting.md#ctrlppexpected-implements-part-of-the-stdexpected-interface) -- it lacks `swap`, `emplace`, the monadic operations and equality, among others.

## Measured figures

- [Real-time safety matrix: estimator stack provenance](rt-safety-matrix.md#provenance) -- the estimator stack figures describe the revision they name and are not re-measured on each later change.
- [Real-time safety matrix: predictive-controller stack provenance](rt-safety-matrix.md#provenance-1) -- the same, for the predictive controllers.
- [Real-time safety matrix: stack cost of the Riccati solve](rt-safety-matrix.md#stack-cost-of-the-riccati-solve) -- the Riccati tables describe the revisions they name.
- [Benchmarks: environment](benchmarks.md#environment) -- the benchmark figures describe the revision they name.
- [Validation: what this page does and does not claim](validation.md#what-this-page-does-and-does-not-claim) -- the cross-validation table records one run in one unpinned environment.
