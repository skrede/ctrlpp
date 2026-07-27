% batch_arx_ident.m -- Batch ARX identification via Octave Control Toolbox arx().
%
% Usage: octave --no-gui batch_arx_ident.m > batch_arx_ident_octave.csv

pkg load control;

% True system: y(k) = 0.8*y(k-1) - 0.2*y(k-2) + 0.5*u(k-1) + 0.3*u(k-2)
a_true = [0.8, -0.2];
b_true = [0.5, 0.3];

na = 2;
nb = 2;
n_samples = 200;

% Generate input-output data
u = zeros(n_samples, 1);
y = zeros(n_samples, 1);

for k = 1:n_samples
    u(k) = sin(0.3 * k) + 0.5 * cos(0.7 * k);
end

for k = 3:n_samples
    y(k) = a_true(1)*y(k-1) + a_true(2)*y(k-2) + b_true(1)*u(k-1) + b_true(2)*u(k-2);
end

% Use arx() from control toolbox
dat = iddata(y, u, 1);
[sys_id, x0] = arx(dat, 'NA', na, 'NB', nb);

% sys_id is a TF model; extract numerator/denominator
[num_id, den_id] = tfdata(sys_id, 'v');

% Some control package versions return cell arrays even for the 'v' (vector) form.
% Unwrap the single-input single-output entry when that happens.
if (iscell(num_id))
    num_id = num_id{1};
endif
if (iscell(den_id))
    den_id = den_id{1};
endif

% ARX model: A(q) y = B(q) u
% den = [1, -a1, -a2, ...], num = [0, b1, b2, ...]
a1 = -den_id(2);
a2 = -den_id(3);
b1 = num_id(2);
b2 = num_id(3);

fprintf('a1,a2,b1,b2\n');
fprintf('%.15e,%.15e,%.15e,%.15e\n', a1, a2, b1, b2);
