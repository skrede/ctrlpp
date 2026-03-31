% lqr_finite_horizon.m -- Finite-horizon LQR via backward Riccati recursion.
% Computes time-varying gains and simulates closed-loop.
%
% Usage: octave --no-gui lqr_finite_horizon.m > lqr_finite_horizon_octave.csv

pkg load control;

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = eye(2);
D_c = zeros(2, 1);
dt = 0.05;

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, dt);
Ad = sys_d.A;
Bd = sys_d.B;

Q = [10 0; 0 1];
R = [1];
Qf = [20 0; 0 2];
horizon = 50;

% Backward Riccati recursion
P = Qf;
gains = zeros(1, 2, horizon);  % K(k) is 1x2

for i = horizon:-1:1
    BtP = Bd' * P;
    S = R + BtP * Bd;
    K = S \ (BtP * Ad);
    gains(:, :, i) = K;
    P = Q + Ad' * P * Ad - Ad' * P * Bd * K;
end

% Simulate closed-loop
x = [1; 0];

fprintf('step,x0,x1,u,K0,K1\n');

for k = 1:horizon
    K = gains(:, :, k);
    u = -K * x;
    fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', k, x(1), x(2), u, K(1), K(2));
    x = Ad * x + Bd * u;
end
