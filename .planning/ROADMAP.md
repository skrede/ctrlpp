# Roadmap: ctrlpp

## Overview

ctrlpp delivers a header-only C++23 control library built on Eigen. PID (with compile-time policy composition), LQR, observers shipped in v0.1.0. Interactive visualization (ctrlpp-implot) shipped in v0.1.1. MPC (linear + nonlinear), comprehensive testing, and composition guide examples shipped in v0.2.0. v0.3.0 shipped nonlinear estimators (EKF, UKF, particle filter), online and offline system identification (RLS, recursive ARX, batch ARX, N4SID), and MPC improvements (general nonlinear constraints, terminal constraint sets). v0.3.1 shipped attitude-aware estimation (MEKF, manifold UKF, complementary filter), signal processing primitives (biquad IIR, FIR, dirty derivative), and moving horizon estimation. v0.3.2 hardens the codebase: removes implot, reorganizes headers into categorical subfolders, adds build tooling, academic references, comprehensive tests, refactoring for readability, gnuplot verification, and world-class hand-written documentation.

## Milestones

- v0.1.0 Foundation, PID, Observers & LQR -- Phases 1-3.1 (shipped 2026-03-15)
- v0.1.1 Interactive Visualization -- Phases 4-9 (shipped 2026-03-18)
- v0.2.0 MPC, Testing & Composition Guide -- Phases 10-13 (shipped 2026-03-21)
- v0.3.0 Nonlinear Estimation, System Identification & MPC Improvements -- Phases 14-17 (shipped 2026-03-22)
- v0.3.1 Estimation & Signal Processing -- Phases 18-21 (shipped 2026-03-23)
- **v0.3.2 Quality, Refactoring & Documentation** -- Phases 22-29 (in progress)

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

<details>
<summary>v0.1.0 Foundation, PID, Observers & LQR (Phases 1-3.1) -- SHIPPED 2026-03-15</summary>

- [x] Phase 1: Foundation & Core Architecture (5/5 plans) -- completed 2026-03-14
- [x] Phase 2: PID Controller (5/5 plans) -- completed 2026-03-15
- [x] Phase 3: Observers & LQR (5/5 plans) -- completed 2026-03-15
- [x] Phase 3.1: Tech Debt Polish (2/2 plans) -- completed 2026-03-15

See [milestones/v0.1.0-ROADMAP.md](milestones/v0.1.0-ROADMAP.md) for full details.

</details>

<details>
<summary>v0.1.1 Interactive Visualization (Phases 4-9) -- SHIPPED 2026-03-18</summary>

- [x] Phase 4: API Refinement (2/2 plans) -- completed 2026-03-16
- [x] Phase 5: Build Infrastructure & Visualization Core (3/3 plans) -- completed 2026-03-16
- [x] Phase 6: Interactive Examples & Static Output (3/3 plans) -- completed 2026-03-17
- [x] Phase 7: Phase 5 Retroactive Verification & CI Fixes (2/2 plans) -- completed 2026-03-17
- [x] Phase 8: Build Path Fix & Pristine Cleanup (3/3 plans) -- completed 2026-03-17
- [x] Phase 9: Runtime & CI Integration Fixes (2/2 plans) -- completed 2026-03-17

See [milestones/v0.1.1-ROADMAP.md](milestones/v0.1.1-ROADMAP.md) for full details.

</details>

<details>
<summary>v0.2.0 MPC, Testing & Composition Guide (Phases 10-13) -- SHIPPED 2026-03-21</summary>

- [x] Phase 10: Linear MPC (5/5 plans) -- completed 2026-03-20
- [x] Phase 11: Nonlinear MPC + NLopt Backend (3/3 plans) -- completed 2026-03-20
- [x] Phase 12: Testing & Hardening (5/5 plans) -- completed 2026-03-20
- [x] Phase 13: Composition Guide Examples (3/3 plans) -- completed 2026-03-21

See [milestones/v0.2.0-ROADMAP.md](milestones/v0.2.0-ROADMAP.md) for full details.

</details>

<details>
<summary>v0.3.0 Nonlinear Estimation, System Identification & MPC Improvements (Phases 14-17) -- SHIPPED 2026-03-22</summary>

- [x] Phase 14: Shared Infrastructure & EKF (3/3 plans) -- completed 2026-03-21
- [x] Phase 15: UKF & Particle Filter (4/4 plans) -- completed 2026-03-21
- [x] Phase 16: MPC Improvements (2/2 plans) -- completed 2026-03-22
- [x] Phase 17: System Identification (3/3 plans) -- completed 2026-03-22

See [milestones/v0.3.0-ROADMAP.md](milestones/v0.3.0-ROADMAP.md) for full details.

</details>

<details>
<summary>v0.3.1 Estimation & Signal Processing (Phases 18-21) -- SHIPPED 2026-03-23</summary>

- [x] Phase 18: SO(3) Primitives, Attitude Estimators & Complementary Filter (4/4 plans) -- completed 2026-03-22
- [x] Phase 19: Signal Processing Primitives (2/2 plans) -- completed 2026-03-22
- [x] Phase 20: Moving Horizon Estimation (2/2 plans) -- completed 2026-03-22
- [x] Phase 21: Testing, Examples & Integration (4/4 plans) -- completed 2026-03-23

See [milestones/v0.3.1-ROADMAP.md](milestones/v0.3.1-ROADMAP.md) for full details.

</details>

### v0.3.2 Quality, Refactoring & Documentation (In Progress)

- [x] **Phase 22: Cleanup & Reorganization** - Remove implot, reorganize flat headers into categorical subfolders (model/, control/, estimation/, lie/), extract shared concepts from mpc/ to model/ (completed 2026-03-24)
- [x] **Phase 23: Tooling & References** - .clang-format, .clang-tidy, CMakePresets.json, coverage option, Doxygen-style citations on all algorithms, REFERENCES.bib (completed 2026-03-24)
- [ ] **Phase 24: Comprehensive Testing** - New unit tests for all types, split oversized files, property-based SO(3) tests, coverage measurement
- [ ] **Phase 25: Refactoring & Gnuplot** - Decompose predict/update into named sub-steps, extract shared helpers, verify gnuplot piping across all examples
- [ ] **Phase 26: Documentation** - Three-tier hand-written docs (getting started, guides, API reference) following mdnspp conventions

## Phase Details

### Phase 22: Cleanup & Reorganization
**Goal**: Codebase contains only live code with categorical header paths -- implot removed, shared concepts in model/, every header in its categorical subfolder
**Depends on**: Nothing
**Requirements**: CLEAN-01, CLEAN-02, CLEAN-03, REORG-01, REORG-02, REORG-03, REORG-04, REORG-05, REORG-06, REORG-07, REORG-08, REORG-09
**Success Criteria** (what must be TRUE):
  1. `lib/ctrlpp-implot/` directory does not exist; no interactive examples remain; `grep -ri implot` returns zero hits (excluding git history and .planning/)
  2. Shared concepts (dynamics_model, measurement_model, differentiable variants) are included from `ctrlpp/model/` by all consumers
  3. Every public header lives in exactly one categorical subfolder (model/, control/, estimation/, lie/, sysid/, mpc/, mhe/, dsp/, detail/)
  4. No forwarding headers exist at old flat paths -- all internal includes use new paths directly; flat convenience headers exist for external users per D-12
  5. Full clean build succeeds; all existing tests pass with new include paths
**Plans:** 3/3 plans complete

Plans:
- [x] 22-01-PLAN.md -- Remove implot library and clean CMake wiring
- [x] 22-02-PLAN.md -- Move headers to categorical subfolders, merge discretise, extract shared concepts, update internal includes
- [x] 22-03-PLAN.md -- Create convenience/umbrella headers, update test/example includes, verify full build

### Phase 23: Tooling & References
**Goal**: Developer tooling configured for the final layout and every algorithm cites its mathematical source
**Depends on**: Phase 22
**Requirements**: BILD-01, BILD-02, BILD-03, REF-01, REF-02, REF-03, REF-04
**Success Criteria** (what must be TRUE):
  1. `clang-format` with the project config produces no changes on the current codebase (or a single formatting pass is applied)
  2. `clang-tidy` runs without errors and enforces snake_case naming on all public identifiers
  3. `cmake --preset dev`, `cmake --preset dev-full`, and `cmake --preset coverage` all configure successfully
  4. Every class implementing a named algorithm has a Doxygen-style `@cite` comment; every citable function cites author, year, equation/section
  5. `REFERENCES.bib` at project root contains 30+ entries covering all algorithm families
**Plans:** 3/3 plans complete

Plans:
- [x] 23-01-PLAN.md -- Create .clang-format, .clang-tidy, and CMakePresets.json tooling configs
- [x] 23-02-PLAN.md -- Create REFERENCES.bib with 30+ academic bibliography entries
- [x] 23-03-PLAN.md -- Add @cite annotations to all algorithm headers and run formatting pass

### Phase 24: Comprehensive Testing
**Goal**: Every public type has focused, human-readable unit tests; oversized test files are split; property-based tests cover SO(3) group axioms; coverage measurement identifies remaining gaps
**Depends on**: Phase 23
**Requirements**: TEST-01, TEST-02, TEST-03, TEST-04, TEST-05, TEST-06, TEST-07, TEST-08, TEST-09, TEST-10, TEST-11, TEST-12, TEST-13
**Success Criteria** (what must be TRUE):
  1. Unit test files exist for SO(3), MEKF, manifold UKF, complementary filter, biquad, FIR, MHE, and NMHE -- each covering edge cases and core functionality
  2. `pid_test.cpp` is split into 3+ focused files; no test file exceeds ~300 lines
  3. Property-based tests verify SO(3) exp/log round-trip, compose associativity, and normalize idempotence
  4. `cmake --preset coverage && cmake --build --preset coverage && ctest --preset coverage` produces an lcov report
  5. All test names are descriptive and readable at a glance (no abbreviations or numeric suffixes)
**Plans:** 4/5 plans executed

Plans:
- [x] 24-01-PLAN.md -- Split pid_test.cpp into 11 focused files by feature cluster
- [x] 24-02-PLAN.md -- Split particle_filter_test, ukf_test, ekf_test into focused files
- [x] 24-03-PLAN.md -- New unit tests for SO(3), biquad, FIR, complementary filter
- [ ] 24-04-PLAN.md -- New unit tests for MEKF, manifold UKF, MHE, NMHE + property-based tests
- [x] 24-05-PLAN.md -- Coverage tooling (codecov.yml) and pipeline verification

### Phase 25: Refactoring & Gnuplot
**Goal**: predict() and update() read like prose, shared logic extracted, and every console example pipes to gnuplot
**Depends on**: Phase 24
**Requirements**: REFAC-01, REFAC-02, REFAC-03, REFAC-04, REFAC-05, REFAC-06, REFAC-07, CLEAN-04, CLEAN-05
**Success Criteria** (what must be TRUE):
  1. predict() methods decomposed into named sub-steps (e.g., propagate_sigma_points, compute_predicted_mean, compute_predicted_covariance)
  2. update() methods decomposed into named sub-steps (e.g., compute_innovation, compute_kalman_gain, apply_state_correction)
  3. All extracted helper functions are `inline`, function templates, or `constexpr` (ODR-safe for header-only)
  4. All existing tests pass with identical numerical results (no tolerance regressions)
  5. Every console example produces gnuplot-pipeable CSV on stdout and has a gnuplot command comment at the top
**Plans:** 3 plans

Plans:
- [ ] 25-01-PLAN.md -- [To be planned]
- [ ] 25-02-PLAN.md -- [To be planned]
- [ ] 25-03-PLAN.md -- [To be planned]

### Phase 26: Documentation
**Goal**: Users can learn, integrate, and reference ctrlpp entirely from hand-written markdown docs -- three tiers from getting started through API reference
**Depends on**: Phase 23, Phase 25
**Requirements**: DOC-01, DOC-02, DOC-03, DOC-04, DOC-05, DOC-06, DOC-07, DOC-08, DOC-09, DOC-10, DOC-11, DOC-12, DOC-13, DOC-14
**Success Criteria** (what must be TRUE):
  1. `docs/README.md` exists with a master index providing three-tier navigation (getting started, guides, API reference)
  2. Getting started guide walks a new user from install through first PID, first observer, and first MPC
  3. Topic guides exist for PID composition, solver injection, observer-controller composition, and sysid workflow
  4. API reference pages exist for every controller, estimator, sysid, DSP, and SO(3) type (~40-50 pages) with working C++23 code examples
  5. Config struct pages document every field with type, default, and description
**Plans:** 3 plans

Plans:
- [ ] 26-01-PLAN.md -- [To be planned]
- [ ] 26-02-PLAN.md -- [To be planned]
- [ ] 26-03-PLAN.md -- [To be planned]

## Progress Table

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation & Core Architecture | v0.1.0 | 5/5 | Complete | 2026-03-14 |
| 2. PID Controller | v0.1.0 | 5/5 | Complete | 2026-03-15 |
| 3. Observers & LQR | v0.1.0 | 5/5 | Complete | 2026-03-15 |
| 3.1 Tech Debt Polish | v0.1.0 | 2/2 | Complete | 2026-03-15 |
| 4. API Refinement | v0.1.1 | 2/2 | Complete | 2026-03-16 |
| 5. Build Infrastructure & Visualization Core | v0.1.1 | 3/3 | Complete | 2026-03-16 |
| 6. Interactive Examples & Static Output | v0.1.1 | 3/3 | Complete | 2026-03-17 |
| 7. Phase 5 Retroactive Verification | v0.1.1 | 2/2 | Complete | 2026-03-17 |
| 8. Build Path Fix & Pristine Cleanup | v0.1.1 | 3/3 | Complete | 2026-03-17 |
| 9. Runtime & CI Integration Fixes | v0.1.1 | 2/2 | Complete | 2026-03-17 |
| 10. Linear MPC | v0.2.0 | 5/5 | Complete | 2026-03-20 |
| 11. Nonlinear MPC + NLopt Backend | v0.2.0 | 3/3 | Complete | 2026-03-20 |
| 12. Testing & Hardening | v0.2.0 | 5/5 | Complete | 2026-03-20 |
| 13. Composition Guide Examples | v0.2.0 | 3/3 | Complete | 2026-03-21 |
| 14. Shared Infrastructure & EKF | v0.3.0 | 3/3 | Complete | 2026-03-21 |
| 15. UKF & Particle Filter | v0.3.0 | 4/4 | Complete | 2026-03-21 |
| 16. MPC Improvements | v0.3.0 | 2/2 | Complete | 2026-03-22 |
| 17. System Identification | v0.3.0 | 3/3 | Complete | 2026-03-22 |
| 18. SO(3), Attitude Estimators & CF | v0.3.1 | 4/4 | Complete | 2026-03-22 |
| 19. Signal Processing Primitives | v0.3.1 | 2/2 | Complete | 2026-03-22 |
| 20. Moving Horizon Estimation | v0.3.1 | 2/2 | Complete | 2026-03-22 |
| 21. Testing, Examples & Integration | v0.3.1 | 4/4 | Complete | 2026-03-23 |
| 22. Cleanup & Reorganization | v0.3.2 | 3/3 | Complete    | 2026-03-24 |
| 23. Tooling & References | v0.3.2 | 3/3 | Complete    | 2026-03-24 |
| 24. Comprehensive Testing | v0.3.2 | 4/5 | In Progress|  |
| 25. Refactoring & Gnuplot | v0.3.2 | 0/TBD | Not started | - |
| 26. Documentation | v0.3.2 | 0/TBD | Not started | - |
