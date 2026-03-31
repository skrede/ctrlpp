% n4sid_ident.m -- N4SID subspace identification via Octave Control Toolbox.
%
% Usage: octave --no-gui n4sid_ident.m > n4sid_ident_octave.csv

pkg load control;

% Generate data from known 2nd-order system
A_true = [0.8, 0.1; -0.2, 0.9];
B_true = [0.5; 0.3];
C_true = [1, 0];
D_true = [0];

n_samples = 500;

% Deterministic excitation
u = zeros(n_samples, 1);
y = zeros(n_samples, 1);
x = [0; 0];

for k = 1:n_samples
    u(k) = sin(0.3 * k) + 0.5 * cos(0.7 * k) + 0.3 * sin(1.1 * k);
    y(k) = C_true * x + D_true * u(k);
    x = A_true * x + B_true * u(k);
end

% Identify using n4sid from control toolbox
dat = iddata(y, u, 1);
sys_id = n4sid(dat, 2);

A_id = sys_id.A;
B_id = sys_id.B;
C_id = sys_id.C;
D_id = sys_id.D;

% Compare via simulation: both models should produce same output
y_pred = zeros(n_samples, 1);
x_id = [0; 0];
for k = 1:n_samples
    y_pred(k) = C_id * x_id + D_id * u(k);
    x_id = A_id * x_id + B_id * u(k);
end

% Output predicted signal (state-space realisation may differ, but output must match)
fprintf('step,y_actual,y_predicted\n');
for k = 1:n_samples
    fprintf('%.15e,%.15e,%.15e\n', k, y(k), y_pred(k));
end
