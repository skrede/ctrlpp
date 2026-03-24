# Requirements: ctrlpp v0.3.2

**Defined:** 2026-03-23
**Core Value:** Lean, composable, and efficient -- users get exactly the controller they need without paying for features they don't use

## v0.3.2 Requirements

Requirements for quality, refactoring & documentation milestone.

### Cleanup

- [x] **CLEAN-01**: Remove lib/ctrlpp-implot entirely (library source, CMake targets, FetchContent deps)
- [x] **CLEAN-02**: Remove all interactive examples that depend on implot
- [x] **CLEAN-03**: Remove implot references from CMake config template, find_package components, and documentation
- [ ] **CLEAN-04**: All console examples produce gnuplot-pipeable CSV output (verified working with gnuplot -p)
- [ ] **CLEAN-05**: Each console example includes a gnuplot command comment at the top showing how to pipe and plot

### Header Reorganization

- [x] **REORG-01**: Shared concepts (dynamics_model, measurement_model, differentiable variants) extracted from mpc/ to core/
- [x] **REORG-02**: Control headers organized into control/ subfolder (pid, lqr, dare, mpc/)
- [x] **REORG-03**: Estimation headers organized into estimation/ subfolder (kalman, luenberger, ekf, ukf, particle_filter, mekf, manifold_ukf, complementary_filter, mhe/)
- [x] **REORG-04**: DSP headers organized into dsp/ subfolder (biquad, fir, cascaded_biquad)
- [x] **REORG-05**: System identification headers remain in sysid/ subfolder (already organized)
- [x] **REORG-06**: Analysis/utility headers organized into core/ (types, state_space, discretize, transfer_function, analysis, so3, observer_policy)
- [x] **REORG-07**: All internal includes updated to new paths -- no forwarding headers (pre-1.0 breaking change is acceptable)
- [x] **REORG-08**: All tests, examples, and CMake updated to new include paths
- [x] **REORG-09**: Build and all tests pass after reorganization

### Testing

- [x] **TEST-01**: Unit tests for SO(3) exp/log/skew/compose/normalize with edge cases (identity, pi rotation, small angle)
- [ ] **TEST-02**: Unit tests for MEKF predict, update, reset, covariance, innovation
- [ ] **TEST-03**: Unit tests for manifold UKF predict, update, geodesic mean, covariance
- [x] **TEST-04**: Unit tests for complementary filter IMU and MARG modes, bias estimation
- [x] **TEST-05**: Unit tests for biquad IIR (low-pass, notch, dirty derivative), frequency response verification
- [x] **TEST-06**: Unit tests for FIR filter, impulse response, convolution correctness
- [ ] **TEST-07**: Unit tests for MHE predict, update, window management, arrival cost
- [ ] **TEST-08**: Unit tests for NMHE predict, update, nonlinear constraints, soft slack
- [x] **TEST-09**: Split pid_test.cpp into focused files (one per policy/feature area)
- [x] **TEST-10**: Split any other test file exceeding ~300 lines into focused files
- [ ] **TEST-11**: Property-based tests for SO(3) (exp/log roundtrip, group axioms, Jacobian consistency)
- [x] **TEST-12**: Coverage measurement integrated via CMake option (CTRLPP_ENABLE_COVERAGE)
- [x] **TEST-13**: All test files are small, focused, with descriptive test names readable at a glance

### Refactoring

- [ ] **REFAC-01**: predict() methods decomposed into named sub-steps (propagate sigma/particles, compute mean, compute covariance, add process noise)
- [ ] **REFAC-02**: update() methods decomposed into named sub-steps (generate measurement predictions, compute innovation, compute cross-covariance, compute Kalman gain, apply correction)
- [ ] **REFAC-03**: Shared helper functions extracted to headers where code is duplicated across estimators
- [ ] **REFAC-04**: Large conditional branches have each body as a separate function
- [ ] **REFAC-05**: Each function does one thing -- no multi-responsibility functions
- [ ] **REFAC-06**: All extracted helpers are inline, template, or constexpr (ODR-safe for header-only)
- [ ] **REFAC-07**: All existing tests pass after refactoring (numerical equivalence verified)

### Academic References

- [x] **REF-01**: Every function implementing citable mathematics has a Doxygen-style comment citing source (author, year, equation/section)
- [x] **REF-02**: Every class implementing a named algorithm has a Doxygen-style header comment with primary reference
- [x] **REF-03**: REFERENCES.bib file at project root with all BibTeX entries (~30-40 expected)
- [x] **REF-04**: References cover: PID (Astrom/Hagglund), LQR/DARE (Anderson/Moore), Kalman (Kalman 1960, Simon 2006), EKF/UKF (Julier/Uhlmann, Wan/van der Merwe), MEKF (Markley/Crassidis), SO(3) (Sola 2018), MPC (Rawlings/Mayne), MHE (Rao 2003), PF (Gordon 1993, Arulampalam 2002), sysid (Ljung, Van Overschee/De Moor), DSP (Oppenheim/Willsky)

### Documentation

- [ ] **DOC-01**: docs/README.md master index with three-tier navigation (getting started -> guides -> API reference)
- [ ] **DOC-02**: Getting started guide (install, first PID, first observer, first MPC)
- [ ] **DOC-03**: Guide: policy-based PID composition
- [ ] **DOC-04**: Guide: concept-based solver injection (MPC/NMPC/MHE)
- [ ] **DOC-05**: Guide: observer-controller composition patterns
- [ ] **DOC-06**: Guide: system identification workflow (identify -> model -> control)
- [ ] **DOC-07**: API reference page for each controller family (PID, LQR, MPC, NMPC)
- [ ] **DOC-08**: API reference page for each estimator (Kalman, Luenberger, EKF, UKF, PF, MEKF, manifold UKF, complementary filter, MHE, NMHE)
- [ ] **DOC-09**: API reference page for sysid (RLS, ARX, N4SID)
- [ ] **DOC-10**: API reference page for DSP (biquad, FIR, cascaded biquad)
- [ ] **DOC-11**: API reference page for SO(3) utilities
- [ ] **DOC-12**: Config struct pages with field tables (type, default, description) for all config types
- [ ] **DOC-13**: All API pages include working C++23 code examples
- [ ] **DOC-14**: Cross-linking between guides, API pages, and related types

### Build & Tooling

- [x] **BILD-01**: .clang-format config matching project style (SortIncludes: Never to respect custom include ordering)
- [x] **BILD-02**: Enhanced .clang-tidy config with readability-identifier-naming for snake_case enforcement
- [x] **BILD-03**: CMakePresets.json for dev, dev-full, and coverage configurations

## Future Requirements

### Ecosystem

- **ECO-01**: vcpkg/Conan packaging recipe
- **ECO-02**: Python bindings (pybind11)
- **ECO-03**: ROS2 wrapper package (ctrlpp_ros)

### Features

- **FEAT-01**: Trajectory/reference generation
- **FEAT-02**: Adaptive control (MRAC/L1)
- **FEAT-03**: SE(3) / Lie group generalization (separate liepp library)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Doxygen HTML generation | Hand-written markdown per mdnspp convention; Doxygen comments in headers for editor tooltips only |
| Static site generator (mkdocs/sphinx) | Pure GitHub markdown rendering is sufficient and zero-maintenance |
| Forwarding headers for old paths | Pre-1.0, breaking changes acceptable; forwarding headers add maintenance burden |
| Benchmark suite | Deferred until feature-complete (per PROJECT.md) |
| CI/CD pipeline changes | Existing CI sufficient for this milestone |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| CLEAN-01 | Phase 22 | Complete |
| CLEAN-02 | Phase 22 | Complete |
| CLEAN-03 | Phase 22 | Complete |
| CLEAN-04 | Phase 25 | Pending |
| CLEAN-05 | Phase 25 | Pending |
| REORG-01 | Phase 22 | Complete |
| REORG-02 | Phase 22 | Complete |
| REORG-03 | Phase 22 | Complete |
| REORG-04 | Phase 22 | Complete |
| REORG-05 | Phase 22 | Complete |
| REORG-06 | Phase 22 | Complete |
| REORG-07 | Phase 22 | Complete |
| REORG-08 | Phase 22 | Complete |
| REORG-09 | Phase 22 | Complete |
| BILD-01 | Phase 23 | Complete |
| BILD-02 | Phase 23 | Complete |
| BILD-03 | Phase 23 | Complete |
| REF-01 | Phase 23 | Complete |
| REF-02 | Phase 23 | Complete |
| REF-03 | Phase 23 | Complete |
| REF-04 | Phase 23 | Complete |
| TEST-01 | Phase 24 | Complete |
| TEST-02 | Phase 24 | Pending |
| TEST-03 | Phase 24 | Pending |
| TEST-04 | Phase 24 | Complete |
| TEST-05 | Phase 24 | Complete |
| TEST-06 | Phase 24 | Complete |
| TEST-07 | Phase 24 | Pending |
| TEST-08 | Phase 24 | Pending |
| TEST-09 | Phase 24 | Complete |
| TEST-10 | Phase 24 | Complete |
| TEST-11 | Phase 24 | Pending |
| TEST-12 | Phase 24 | Complete |
| TEST-13 | Phase 24 | Complete |
| REFAC-01 | Phase 25 | Pending |
| REFAC-02 | Phase 25 | Pending |
| REFAC-03 | Phase 25 | Pending |
| REFAC-04 | Phase 25 | Pending |
| REFAC-05 | Phase 25 | Pending |
| REFAC-06 | Phase 25 | Pending |
| REFAC-07 | Phase 25 | Pending |
| DOC-01 | Phase 26 | Pending |
| DOC-02 | Phase 26 | Pending |
| DOC-03 | Phase 26 | Pending |
| DOC-04 | Phase 26 | Pending |
| DOC-05 | Phase 26 | Pending |
| DOC-06 | Phase 26 | Pending |
| DOC-07 | Phase 26 | Pending |
| DOC-08 | Phase 26 | Pending |
| DOC-09 | Phase 26 | Pending |
| DOC-10 | Phase 26 | Pending |
| DOC-11 | Phase 26 | Pending |
| DOC-12 | Phase 26 | Pending |
| DOC-13 | Phase 26 | Pending |
| DOC-14 | Phase 26 | Pending |

**Coverage:**
- v0.3.2 requirements: 55 total
- Mapped to phases: 55
- Unmapped: 0

---
*Requirements defined: 2026-03-23*
*Last updated: 2026-03-23 after roadmap creation*
