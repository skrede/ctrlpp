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

`ctrlpp::expected` carries `[[nodiscard]]` at class level (see the two class
declarations in `ctrlpp/detail/expected.h`), on both the
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

These queries carry **no** `[[nodiscard]]`. Discarding one is
harmless: it is a
question the caller may ask whenever it wants, and never asking it is a legitimate
choice.

## 4. Where the attribute goes, and where it does not

`[[nodiscard]]` is written at class level on the result type and on
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
Every solver backend adapter reports setup failure through a fallible
`setup(problem)` returning `ctrlpp::expected<void, E>`, and the dispatchers
`setup_qp_solver` and `setup_nlp_solver` forward that typed error rather than
flattening it, so the cause of a setup failure is readable at the seam. The
moving-horizon and nonlinear moving-horizon estimators report their EKF fallback
through `diagnostics()` (channel 2).

The per-step measurement `update` on all seven estimator types is on channel 1:
`kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf`, `luenberger_observer`,
and all three `complementary_filter` overloads return
`ctrlpp::expected<void, E>` with a per-module error enumeration
(`kalman_update_error`, `ekf_update_error`, and so on -- one per module, not one
shared enum). Each rejects a non-finite operand **before mutating anything**, so
a rejected step leaves the carried estimate bitwise unchanged. Each also carries
the channel-3 state-health query: `health()`, returning a per-module latching
enumeration, so a caller can ask whether the estimate it is carrying is still
degraded from an earlier step. Two of those queries predate this convention
(`ukf_health`, `manifold_ukf_health`) and were extended rather than replaced.
None of the seven carries a discard annotation on `health()`.

The complementary filter's zero-norm acceleration and magnetic skips stay on the
success path and are **not** rejections. They are channel-2 shaped -- the step
ran and produced a valid attitude, it simply had no correction to apply -- and
they are candidates for an explicit disposition report rather than for channel 1.

The per-step `compute` and `evaluate` surfaces on the three stateful controllers
are on channel 1 as well: `pid::compute` (both overloads),
`mrac_controller::evaluate` and `l1_controller::evaluate` return
`ctrlpp::expected<vector_t, E>` -- the result **carries the control vector**,
because a controller step produces the command it was asked for and channel 1 is
the result type over whatever the operation produces. Each carries a per-module
enumeration (`pid_step_error`, `mrac_step_error`, `l1_step_error`), rejects
before mutating anything, and carries an unannotated channel-3 `health()` query
(`pid_health`, `mrac_health`, `l1_health`).

A refusal on these three is not the same event as a refusal on an estimator
update. A refused estimator step leaves the estimate standing; a refused
controller cycle means the caller has **no control output for this cycle**, and
an actuator is going to be driven by something regardless. The caller chooses
what -- hold the last successful command, drive a configured safe value, or fail
over -- and the API deliberately does not choose, because the right answer is a
property of the plant. This is also why none of the three was reduced to
`expected<void, E>` with the command written through a reference parameter: that
would satisfy the shape of the convention while deleting the operation's output.

`l1_controller` additionally reports on channel 3 a success-path fact that no
per-cycle result could carry: its clamping projection substitutes a configured
bound for an infinite adaptation estimate, producing a finite in-range command
from a meaningless estimate. `health()` latches
`projection_clamped_non_finite` when that happens.

### Configuration validated at construction

`kalman_filter`, `ekf`, `ukf` and `rls` now validate their configuration on
channel 1 at construction, so a mistake made before the object ever ran is
reported where it was made rather than as a non-finite estimate at the first
step. Two of the four had **no factory at all** and gained one; a third had a
`try_create`; only `rls` had an infallible construction to convert.

Each of the four exposes `create` as its **only** construction path -- the plain
constructors are private, reachable solely through a `validated_tag` the factory
holds. A validating factory beside a public constructor validates nothing,
because any caller can take the other door. That also settled the naming: the
`try_` prefix marks a fallible factory that contrasts with a genuinely
non-fallible constructor, so with `ukf`'s constructors gone the prefix
contrasted with nothing and the member became `create`. The sigma-point
strategies keep `try_create`, because their own plain constructors do survive.

The three filters share `filter_error`, which gained
`non_finite_process_noise`, `non_finite_measurement_noise`,
`non_finite_initial_state` and `non_finite_initial_covariance` -- one shared
enumeration rather than three identical copies, because the three configuration
aggregates declare the same four fields and feed them into the same recursion.
`rls` carries its own `rls_error`, whose conditions are read off the covariance
update rather than asserted as a range.

**Finiteness is the domain condition and the whole of it.** An ill-conditioned
but finite configuration -- a covariance with entries many orders of magnitude
apart, a singular `P0` -- is accepted. Rejecting it would convert a
numerical-behavior question into a domain violation and refuse problems the
library solves. `rls_error` keeps one enumerator that is explicitly *not* a
domain condition, `forgetting_factor_above_unity`, and says so: the arithmetic
there is well defined and the consequence is an estimator that silently stops
adapting, so it is enforced as the type's stated `(0, 1]` contract, not as
arithmetic.

The two moving-horizon estimators embed an extended filter and hand it the same
noise fields, so `mhe` and `nmhe` became fallible too and forward that filter's
rejection verbatim. That is not plumbing: both also invert `Q` and `R` to form
the arrival-cost and stage weights, so a non-finite entry poisons the posed
problem as well as the filter.

Not yet converted: `predict` on those seven estimator types; the configuration
and reset paths of the three controllers (`pid::set_params`, `pid::set_integral`,
the adaptive controllers' construction and `reset`), through which a non-finite
value still enters silently and is caught only by the next cycle's rejection and
the `health()` latch; `rls::update`, which still returns a bare `bool` and so
cannot say why a sample was skipped; `ekf_config::numerical_eps`, whose own
domain condition (finite and strictly positive, since the central-difference
stencil divides by it) is not yet checked at construction; the configuration
validation of `mekf`, `manifold_ukf` and `complementary_filter`, whose factories
exist but validate only their initial quaternion and not their `Q`, `R` or `P0`;
`particle_filter::update` and the online planners' `update`; and `lqr_gain` and
`lqi_gain`, which still return `std::optional` and so discard the reason for the
empty result. Those surfaces are being moved onto
channel 1. Until each one is, the NaN and Inf propagation contract described in
[numerical-behavior.md](numerical-behavior.md) is what governs them.
