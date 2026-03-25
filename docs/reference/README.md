# Reference

Standalone theory and mathematical background for the algorithms implemented in
ctrlpp. Each page presents key equations in LaTeX notation, explains the
underlying theory, and cites the academic sources from
[REFERENCES.bib](../../REFERENCES.bib).

API pages link here for derivations and context; theory pages link back to their
corresponding API pages.

## Pages

- [PID Theory](pid-theory.md) -- parallel form, derivative filter, anti-windup clamping and back-calculation
- [Kalman Theory](kalman-theory.md) -- linear Kalman filter predict/update equations, optimality, innovation
- [EKF/UKF Theory](ekf-ukf-theory.md) -- EKF Jacobian linearisation, UKF sigma-point generation and weights, unscented transform
- [Particle Filter Theory](particle-filter-theory.md) -- importance sampling, weight update, systematic resampling, ESS
- [Attitude Estimation Theory](attitude-estimation-theory.md) -- quaternion kinematics, MEKF error-state formulation, Rodrigues' formula, manifold UKF
- [MPC Theory](mpc-theory.md) -- QP/NLP optimisation formulation, terminal cost and constraints, stability
- [MHE Theory](mhe-theory.md) -- moving horizon cost function, arrival cost approximation, duality with MPC
- [System Identification Theory](sysid-theory.md) -- ARX regression model, RLS recursive update, N4SID Hankel matrix decomposition
- [DSP Theory](dsp-theory.md) -- biquad transfer function, bilinear transform, FIR convolution, cascading
