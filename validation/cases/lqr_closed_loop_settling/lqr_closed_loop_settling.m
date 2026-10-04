% lqr_closed_loop_settling.m -- The embedded golden cell, answered by the Octave Control Toolbox.
%
% Usage: octave --no-gui lqr_closed_loop_settling.m > lqr_closed_loop_settling_octave.csv

pkg load control;

% Both values restate examples/embedded/shared/golden_reference.h, which the C++ arm
% reads; a drift there turns this case red rather than passing unnoticed.
dt = 0.02;
steps = 201;

A = [1 dt; 0 1];
B = [0.5 * dt^2; dt];
Q = diag([10 1]);
R = 0.1;

[P, ~, ~] = dare(A, B, Q, R);

% Written out rather than taken from a helper: this is the gain for u = -K*x,
% the sign convention the C++ arm's step() applies.
K = (R + B' * P * B) \ (B' * P * A);

x = [1; 0];
for k = 1:steps
  u = -K * x;
  x = A * x + B * u;
end

fprintf('K_00,K_01,x0_final,x1_final,final_norm\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e\n', K(1,1), K(1,2), x(1), x(2), norm(x));
