% dare_solution.m -- Reference DARE solution via Octave Control Toolbox.
% Outputs the upper-triangular elements of P as a single-row CSV (reproducible, no time series).
%
% Usage: octave --no-gui dare_solution.m > dare_solution_octave.csv

pkg load control;

% Mass-spring-damper: m=1, k=1, b=0.5, dt=0.05
% Continuous: A_c = [0 1; -1 -0.5], B_c = [0; 1]
% Discretise with ZOH via c2d, then solve DARE.

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = eye(2);
D_c = zeros(2, 1);

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, 0.05);

A = sys_d.A;
B = sys_d.B;
Q = [10 0; 0 1];
R = [1];

[P, ~, ~] = dare(A, B, Q, R);

% Also compute LQR gain K = (R + B'PB)^{-1} B'PA
K = (R + B' * P * B) \ (B' * P * A);

% Output P as flattened row-major, and K as flattened row-major
% This avoids ambiguity about matrix layout.
fprintf('P_00,P_01,P_10,P_11,K_00,K_01\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        P(1,1), P(1,2), P(2,1), P(2,2), K(1,1), K(1,2));
