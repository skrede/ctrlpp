% care_solution.m -- Reference CARE solution via Octave Control Toolbox.
% Outputs P (row-major) and K (row-major) as a single CSV row.
%
% Usage: octave --no-gui care_solution.m > care_solution_octave.csv

pkg load control;

% Mass-spring-damper: m=1, k=1, b=0.5
% Continuous: A_c = [0 1; -1 -0.5], B_c = [0; 1].
% Solve continuous-time CARE A'P + PA - PBR^{-1}B'P + Q = 0 directly.

A = [0, 1; -1, -0.5];
B = [0; 1];
Q = [10 0; 0 1];
R = [1];

[P, ~, ~] = care(A, B, Q, R);

% Continuous-time LQR gain K = R^{-1} B' P
K = R \ (B' * P);

fprintf('P_00,P_01,P_10,P_11,K_00,K_01\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        P(1,1), P(1,2), P(2,1), P(2,2), K(1,1), K(1,2));
