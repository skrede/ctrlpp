# Error Reporting

ctrlpp reports through two channels, and answers a third kind of question outside
both of them. Which one an operation uses is decided by a single question: did the
operation do what the caller asked?

## 1. Failure

A fallible operation returns `ctrlpp::expected<T, E>`, where `E` is a per-module
error enumeration. This is the channel for "I could not do what you asked":
construction, setup, solve, update, reconfigure.

```cpp
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          detail::conditioning_policy Cond = detail::pivot_ratio_conditioning>
auto dare(const Eigen::Matrix<Scalar, int(NX), int(NX)>& A,
          const Eigen::Matrix<Scalar, int(NX), int(NU)>& B,
          const Eigen::Matrix<Scalar, int(NX), int(NX)>& Q,
          const Eigen::Matrix<Scalar, int(NU), int(NU)>& R,
          Cond                                           /*tag*/ = {})
    -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>;
```

The caller branches on the result and reaches the value through `operator*`:

```cpp
auto ctrl_result = controller::create(cfg, 5.0, 100.0);
if(!ctrl_result.has_value())
{
    std::cerr << "invalid L1 configuration\n";
    return 1;
}
auto& ctrl = *ctrl_result;
```

A rejection is `return ctrlpp::unexpected(<enum>::<value>);`. Success for an
operation that produces no value is `return {};`.

`ctrlpp::expected` carries the standard discard-warning attribute at class level
(see the two class declarations in `ctrlpp/detail/expected.h`), on both the
primary template and the `void` partial specialization, so ignoring a failure is
a compile-time diagnostic rather than a silent runtime no-op. ctrlpp's own build
promotes that diagnostic to an error on its own targets, so a discarded fallible
return does not build here at all. Because the attribute sits on the type, every
fallible return in the library inherits it; nothing has to be annotated per
function, and nothing can drift out of sync.

## 2. Disposition

When an operation **succeeded** but produced something other than what was
commanded, that is not a failure and is not reported through `ctrlpp::expected`.
A fallback taken, a substituted shape, a limit raised to make the request
feasible: the operation ran, and it produced a valid, usable result. The report is
a separate aggregate the caller reads back from the object.

The reason for keeping this out of the failure channel is behavioral, not
stylistic. Forcing a caller to unwrap a non-failure through the failure path
teaches the caller that the failure path is usually noise, and a caller who has
learned that will eventually ignore a real failure.

The moving-horizon estimator is the worked example. Its diagnostics aggregate
lives in `ctrlpp/mhe/mhe_diagnostics.h`:

```cpp
template <typename Scalar>
struct mhe_diagnostics
{
    solve_status status{solve_status::error};
    int iterations{};
    Scalar solve_time{};
    Scalar cost{};
    Scalar primal_residual{};
    Scalar dual_residual{};
    Scalar max_constraint_violation{};
    Scalar max_residual_bound_violation{};
    Scalar total_slack{};
    bool used_ekf_fallback{false};
};
```

and is read back through `diagnostics()` in `ctrlpp/mhe.h`:

```cpp
const mhe_diagnostics<Scalar>& diagnostics() const { return m_diagnostics; }
```

When the estimator falls back to an EKF step, it produced a valid estimate. It
did not fail. `used_ekf_fallback` plus `status` is how the caller learns which of
the two it got.

## 3. State health

A persistent question about the object, of the form "is the state I am carrying
still degraded from something that happened earlier", is neither of the above. It
is a plain `const` member function returning a per-module status enumeration.

A per-call disposition cannot answer it, because the question outlives the call:
a covariance that had to be repaired three steps ago is still a repaired
covariance now, and the call that just returned cleanly has nothing to say about
it.

The unscented filter is the in-tree precedent. From
`ctrlpp/estimation/ukf.h`:

```cpp
enum class ukf_health
{
    ok,
    covariance_repaired
};
```

```cpp
/// @brief Report whether the filter has had to repair a non-positive-definite
/// covariance to the nearest symmetric positive definite matrix.
ukf_health health() const { return m_health; }
```

These queries carry **no** discard-warning attribute. Discarding one is
harmless: it is a
question the caller may ask whenever it wants, and never asking it is a legitimate
choice.

## 4. Where the attribute goes, and where it does not

The discard-warning attribute is written at class level on the result type and on
the disposition aggregates, and nowhere else in the library. It is never written
at a call site or on an individual function declaration.

It is **forbidden** on:

- accessors and getters
- size and empty queries
- state-health predicates
- any other return whose discard is harmless

The attribute is earned where discarding the result silently breaks correctness.
Dropping a reported failure inside a control loop leaves a step that did nothing
while the caller believes it succeeded, and destroys the enumerator that named
the cause; dropping a substitution report leaves the caller acting on a command
the library quietly changed. Those are real misuse, and stopping them at compile
time is worth the ceremony it costs.

Putting the same attribute on a trivial return buys none of that and charges the
caller for it anyway. A caller who wants the size and then decides not to use it
has done nothing wrong, and making them write a cast to say so is noise that
trains them to write the cast reflexively, including at the one call site where it
mattered.

## Current state of the conversion

**This describes the library as it stands at this commit, not a finished state.**

Following the convention today: the fallible `create` factories and the discrete
algebraic Riccati solver report failure through `ctrlpp::expected` (channel 1).
The moving-horizon and nonlinear moving-horizon estimators report their EKF
fallback through `diagnostics()` (channel 2). The unscented filter exposes
`health()` as a state-health query.

Not yet converted: the per-step `update`, `compute`, and `evaluate` surfaces on
the filters and controllers still return `void` and report nothing at all; the
solver setup dispatchers `setup_qp_solver` and `setup_nlp_solver` still return
`bool`; and `lqr_gain` and `lqi_gain` still return `std::optional`, which
discards the reason for the empty result. Those surfaces are being moved onto
channel 1. Until each one is, the NaN and Inf propagation contract described in
[numerical-behavior.md](numerical-behavior.md) is what governs them.
