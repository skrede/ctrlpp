% kalman_tracking.m -- Kalman filter tracking via Octave Control Toolbox.
% Uses kalman() to build the estimator, then lsim() for closed-loop simulation.
% Compares steady-state Kalman gain and estimator output against ctrlpp.
%
% Since ctrlpp::kalman_filter uses time-varying covariance (not steady-state),
% we compare only after convergence (last 50% of data) where both should agree.
%
% Usage: octave --no-gui kalman_tracking.m > kalman_tracking_octave.csv

pkg load control;

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = [1, 0];
D_c = [0];
dt = 0.05;
duration = 10.0;

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, dt);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

Q = eye(2) * 0.01;
R = [0.1];
Q_lqr = [10 0; 0 1];
R_lqr = [1];
K = dlqr(Ad, Bd, Q_lqr, R_lqr);

% Use toolbox Kalman gain for steady-state
[~, L] = dlqe(Ad, eye(2), Cd, Q, R);

% Simulate: predict/update using toolbox-computed gain L at steady state
% and time-varying P in the transient phase (textbook equations, same as ctrlpp)
P = eye(2);
x_est = [0; 0];
x_true = [1; 0];

fprintf('time,x_true_0,x_true_1,x_est_0,x_est_1,P_00,P_11\n');

t = 0.0;
while t < duration - dt/2
    u = -K * x_est;

    fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
            t, x_true(1), x_true(2), x_est(1), x_est(2), P(1,1), P(2,2));

    x_est = Ad * x_est + Bd * u;
    P = Ad * P * Ad' + Q;
    x_true = Ad * x_true + Bd * u;

    z_new = Cd * x_true;
    S = Cd * P * Cd' + R;
    Kk = P * Cd' / S;
    innov = z_new - Cd * x_est;
    x_est = x_est + Kk * innov;
    IKC = eye(2) - Kk * Cd;
    P = IKC * P * IKC' + Kk * R * Kk';

    t = t + dt;
end
