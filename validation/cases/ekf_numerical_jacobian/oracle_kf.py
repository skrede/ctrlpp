#!/usr/bin/env python3
"""Independent scipy/numpy oracle for the linear constant-velocity Kalman
recursion exercised by the "ekf with numerical Jacobians converges on linear
system" test case (tests/unit/ekf_convergence_test.cpp).

The dynamics and measurement in that test case are linear, so the analytic
Jacobians F, G, H are exact closed forms; the EKF's central-difference
numerical Jacobian should reproduce them up to the finite-difference
round-off floor. This script implements the identical linear Kalman
recursion (Joseph-form covariance update, matching
ctrlpp/estimation/ekf.h::update_covariance) using only numpy, independent of
ctrlpp, and prints the state estimate after 50 predict/update steps to full
double precision.

Run: python3 oracle_kf.py
Cross-checked against oracle_kf.m (Octave); see the agreement recorded in
the SUMMARY for this plan.
"""
import numpy as np

dt = 0.1
F = np.array([[1.0, dt], [0.0, 1.0]])
G = np.array([[0.5 * dt * dt], [dt]])
H = np.array([[1.0, 0.0]])
Q = np.eye(2) * 0.01
R = np.array([[1.0]])

x = np.zeros((2, 1))
P = np.eye(2) * 10.0

true_pos = 0.0
true_vel = 1.0

for i in range(50):
    true_pos += true_vel * dt
    u = 0.0

    # predict
    x = F @ x + G * u
    P = F @ P @ F.T + Q

    # measurement
    z = true_pos + 0.1 * np.sin(float(i))

    # update (Joseph form)
    z_pred = (H @ x)[0, 0]
    innovation = z - z_pred
    S = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(S)
    x = x + K * innovation
    IKH = np.eye(2) - K @ H
    P = IKH @ P @ IKH.T + K @ R @ K.T

print(f"final_pos = {x[0, 0]!r}")
print(f"final_vel = {x[1, 0]!r}")
print(f"final_P00 = {P[0, 0]!r}")
print(f"final_P11 = {P[1, 1]!r}")
print(f"true_pos  = {true_pos!r}")
print(f"true_vel  = {true_vel!r}")
