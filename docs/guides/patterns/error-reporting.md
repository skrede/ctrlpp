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
    -> ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error>;
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

A refused measurement is **not** on this channel, and putting it here was the
mistake worth naming. Nothing was produced for the aggregate to describe, and the
one accessor a caller would reach for -- `state()` -- goes on returning the
estimate the last accepted measurement produced. So `mhe::update` and
`nmhe::update` return `expected<void, ekf_update_error>`, forwarding the embedded
filter's verdict, and the aggregate goes on describing the last step that
succeeded. The two accessors then always describe the same step, and
`used_ekf_fallback == false` means one thing -- the window solve produced this
estimate -- rather than doubling as "nothing produced anything".

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
    covariance_repaired,
    non_finite_estimate
};
```

```cpp
/// @brief Report whether the carried estimate is still degraded from an
/// earlier step: a covariance that was not positive definite, or a
/// non-finite estimate. Latches; a rejected measurement does not set it.
ukf_health health() const { return m_health; }
```

The status is latched where the degradation happens, not where it is first
consumed: `ukf::predict` and `manifold_ukf::predict` judge the estimate they
produced and latch the rung it warrants, so a caller reading `covariance()`
between `predict` and `update` never sees an indefinite or non-finite
prediction beside `ok`.

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

`particle_filter::update` is on channel 1 as well, returning
`ctrlpp::expected<void, pf_update_error>` (`non_finite_state`,
`non_finite_measurement`, `non_finite_result`). It refuses before committing a
weight or drawing from its random generator, so a refused step leaves the
particles, the weights and the generator bitwise unchanged. Its channel-3 query
`health()` returns `pf_health` (`ok`, `weights_reset`,
`resampling_rolled_back`, `non_finite_estimate`). A weight depletion, where no
particle explains a finite measurement, is not a rejection: the step resets the
weights to uniform, succeeds with a valid but uninformed posterior, and latches
`weights_reset`, so it stays distinguishable from a refused measurement. A
rolled-back resampling is the same kind of success-path rung: when roughening
would leave a particle non-finite, the step restores the particles from before
the resampling, keeps its conditioned weights, succeeds, and latches
`resampling_rolled_back`.

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

A variant keeps the factory the only door while letting the result be built in
place: the constructor is public, but its tag is a passkey whose default
constructor is private and which only the class can make. `ctrlpp::expected`
can then construct the object in its own storage through its `std::in_place`
constructor, and guaranteed copy elision carries that storage to the caller, so
the factory never holds a whole object in its frame. `particle_filter` uses it,
because its object holds the whole particle set.

The three filters share `filter_error`, which gained
`non_finite_process_noise`, `non_finite_measurement_noise`,
`non_finite_initial_state` and `non_finite_initial_covariance` -- one shared
enumeration rather than three identical copies, because the three configuration
aggregates declare the same four fields and feed them into the same recursion.
`rls` carries its own `rls_error`, whose conditions are read off the covariance
update rather than asserted as a range.

**Every filter that takes a `Q`, `R` or `P0` refuses one that is non-finite,
structurally asymmetric or indefinite.** A matrix with an eigenvalue below
`-floor`, with the floor `max(N eps max|M_ij|, N denorm_min)` (the symmetric
eigensolver's backward-error order, kept from underflowing to zero), is not a
covariance, whatever the filter does with it, so `filter_error` also carries
`indefinite_process_noise`, `indefinite_measurement_noise` and
`indefinite_initial_covariance`. Nor is a matrix whose triangles differ by more
than rounding: `filter_error` carries `asymmetric_process_noise`,
`asymmetric_measurement_noise` and `asymmetric_initial_covariance` for an
asymmetry `max|M_ij - M_ji|` beyond the counted bound
`max(N eps max|M_ij|, N denorm_min)`. An asymmetry within that bound is drift:
the matrix is accepted, the filter stores its symmetrized form (each differing
pair replaced by its midpoint), and the definiteness check runs on that
symmetrized matrix. Each matrix is checked for finiteness, then for symmetry,
then for definiteness, and the matrices in declaration order; the convention
and its count are in
[numerical-behavior.md](numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused).
A **singular** positive
semidefinite matrix is accepted: a point-mass `P0`, a rank-one `Q` and a
noiseless `R` are legitimate models, and the Kalman-family recursions invert
only the innovation covariance, never `R` alone. An ill-conditioned
configuration, such as a covariance with entries many orders of magnitude
apart, is accepted too: conditioning is a numerical-behavior question, not a
domain violation. `rls_error` keeps one enumerator that is explicitly *not* a
domain condition, `forgetting_factor_above_unity`, and says so: the arithmetic
there is well defined and the consequence is an estimator that silently stops
adapting, so it is enforced as the type's stated `(0, 1]` contract, not as
arithmetic.

`mekf` and `manifold_ukf` apply the same checks to their `Q`, `R` and `P0` after
the initial quaternion, and `mekf` refuses a non-finite initial bias between `R`
and `P0`, the order its configuration declares them, with its own
`non_finite_initial_bias`: the bias seeds part of the carried state, but it is
not the initial state `x0` the other filters take. All three attitude factories,
`mekf`, `manifold_ukf` and `complementary_filter`, check the initial quaternion
first: a zero or non-finite norm is `degenerate_quaternion`, and a norm that
differs from one by more than the rounding bound `max(4 eps, 4 denorm_min)` is
`non_unit_quaternion`, because a rotation quaternion is unit by definition. A
norm defect within the bound is drift and is normalized away, under the same
[convention](numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused)
as a covariance's asymmetry. `kalman_filter::set_noise` and
`kalman_filter::reset_covariance`, the two setters that replace a covariance
after construction, return `ctrlpp::expected<void, filter_error>` and refuse
with the same enumerators before assigning anything; `set_noise` validates both
matrices before it assigns either.

The two moving-horizon estimators embed an extended filter and hand it the same
noise fields, so `mhe` and `nmhe` became fallible too and forward that filter's
rejection verbatim, the symmetry and definiteness refusals included. That is not plumbing:
both also invert `Q` and `R` to form the arrival-cost and stage weights, so a
non-finite entry poisons the posed problem as well as the filter, and a
singular `Q` or `R`, which the filter accepts, is refused by their own
`moving_horizon_configuration_error`.

`particle_filter` validates its configuration on channel 1 as well, through
`particle_filter::create` and `make_particle_filter`, which return
`ctrlpp::expected<particle_filter, filter_error>`; its constructor takes the
`validated_tag` passkey only the class can make, and the factory builds the
filter in the caller's storage. It refuses the same finiteness, symmetry and
definiteness conditions, samples and scores the symmetrized matrices, and, in
addition, refuses a singular `R` (`singular_measurement_noise`): it
scores every particle with R's inverse and log-determinant, which a singular R
does not have. A singular but positive semidefinite `Q` or `P0` is still
accepted, because a rank-one process noise and a point-mass prior have a
Gaussian to sample.

`pid` validates its configuration on channel 1: `pid::create`,
`pid::set_params` and `pid::set_integral` refuse a value that cannot work with a
`pid_config_error` and leave the controller untouched. `create` is its only
construction path; the constructor is private.

Not yet converted: `predict` on those seven estimator types (the two unscented
filters' `predict` latches `health()` on what it produced, but still refuses
nothing); the configuration
and reset paths of the two adaptive controllers (their construction and
`reset`), through which a non-finite value still enters silently and is caught
only by the next cycle's rejection and the `health()` latch; `rls::update`, which still returns a bare `bool` and so
cannot say why a sample was skipped; `ekf_config::numerical_eps`, whose own
domain condition (finite and strictly positive, since the central-difference
stencil divides by it) is not yet checked at construction; the remaining
configuration fields of the quaternion filters, `mekf_config::dt` and
`numerical_eps`, the geodesic-mean options of `manifold_ukf_config`, and the
gains and step of `cf_config`, whose `complementary_filter` factory validates
only the initial quaternion; and the online planners' `update`. Those surfaces are
being moved onto channel 1. Until each one is, the NaN and Inf propagation
contract described in [numerical-behavior.md](numerical-behavior.md) is what
governs them.

### No failure is signalled by an empty optional any more

`std::optional` no longer carries a failure anywhere in `lib/`. Eight surfaces
moved onto channel 1:

| Surface | Enumeration | What the empty result used to discard |
| --- | --- | --- |
| `lqr_gain` (both discrete overloads) | `dare_error`, **forwarded** | which of six conditions refused the pair |
| `lqr_gain_continuous` | `care_error`, forwarded | the same, plus its own non-finite rejection |
| `detail::partition_lqi_gain`, `lqi_gain` | `dare_error`, forwarded | the augmented solve's cause |
| `place`, `place_observer` | `place_error`, **new** | a structural refusal versus a numerical one |
| `so3::normalize` | `so3_error`, new | that the postcondition was not met at all |
| `dsp::detail::validate_biquad_design` | `dsp_error` | nothing -- see below |

The three LQR forms and both LQI forms **forward** an enumeration that already
existed one call below rather than inventing a name for it. The cause was
computed and typed by the Riccati solver and then thrown away at the seam; the
repair is to stop throwing it away, not to describe it a second time where the
two descriptions can drift apart.

`place` had **no** enumeration, so one was derived from its actual refusal
sites, one enumerator per distinct cause. The dividing line that mattered there
is structural versus numerical: `multi_input_not_supported` and
`multi_output_not_supported` cannot be lifted by changing any number, while
`poles_not_conjugate_symmetric` and `uncontrollable_pair` are conditions on the
data. Sharing an enumerator across that line would send a caller to change the
one thing that cannot help.

`validate_biquad_design` is the exception that proves the rule and is recorded
as such: it returned `std::optional<dsp_error>`, an **inverted** channel where
an empty result meant success and an engaged one carried the error. It lost no
information. It was converted so the library has one channel shape rather than
two, and that is its only justification: it is not a defect fix.

Every call site branches. None substitutes a default value for a refused
result, because relocating a dishonesty is not removing it.

### Two enumerators added, for causes that were being reported as something else

Both Riccati enumerations gained `singular_r`. Their pencil and Hamiltonian
builds invert the input weighting exactly as they invert the state matrix, and
`singular_a` already existed for the latter; the missing counterpart meant a
caller who set a weighting to zero deliberately was told their input was
non-finite. Worse, a weighting that is **rank-deficient but nonzero** never went
non-finite at all -- the rank-revealing solve returns a least-squares answer over
the leading rank columns -- so both solvers ran to completion and reported
success on a problem the caller had not posed. The condition is read off the
factorization each build already forms, so no second factorization was added.

`singular_u11` in both enumerations was **documented rather than renamed**, which
is a distinction worth keeping straight. It covers two situations at once: a pair
with no stabilizing solution that the eigenvalue-count test cannot see, and a
numerical failure to separate the invariant subspace on a genuinely well-posed
pair. Both were measured on one sweep. Renaming it after the structural cause
would over-claim on every well-posed row; leaving it undocumented would leave a
caller reading a numerical symptom as the whole story. Stating what it covers,
and that it does not distinguish the two, is the honest option and costs nothing
at runtime.

That is not the same as writing a comment in place of a fix. Where a cause is
already computed and typed one call below, forwarding it is free and a comment
explaining the loss would be the defect wearing a fix's clothes. Here **no
mechanism computes the cause at all**, and adding one is a stabilizability test
with its own threshold to derive.

### An honest postcondition is a channel-1 question too

`so3::normalize` documented a unit-norm postcondition and did not deliver it for
three families of input, with no channel to say so. It is now fallible, with
`so3_error::non_finite_input` and `so3_error::zero_quaternion` naming the two
inputs that genuinely have no unit representative; every other finite quaternion
is normalized, including ones whose squared norm is not representable. See
[so3](../../api/lie/so3.md).

### A limit the caller stated is refused by default, and raised only on request

`trapezoidal_trajectory` cannot reconcile some boundary velocities over the
commanded displacement at the commanded acceleration. By default that is a
failure, `trajectory_error::acceleration_limit_infeasible`, because returning a
profile above a limit the caller stated is the silent success channel 1 exists
to prevent. A caller whose limit is advisory opts in with
`acceleration_limit_policy::raise_to_feasible`; the operation then **succeeds**
with the profile correct under the raised limit, and the raise is channel 2,
reported through `disposition()` -- carrying the commanded and the realized
acceleration, not a flag. See
[trapezoidal_trajectory](../../api/trajectory/trapezoidal-trajectory.md).
