% tf_ss_conversion.m -- TF<->SS round-trip via Octave.
% Verifies that ss2tf(tf2ss(H)) recovers original TF coefficients.
% Also compares poles (eigenvalues of A), which are invariant across canonical forms.
%
% Usage: octave --no-gui tf_ss_conversion.m > tf_ss_conversion_octave.csv

pkg load control;
pkg load signal;

% Transfer function: H(s) = (2s + 3) / (s^2 + 0.5s + 1)
num = [2, 3];
den = [1, 0.5, 1];

% TF -> SS
[A, B, C, D] = tf2ss(num, den);

% Poles from A matrix
p = sort(eig(A));

% SS -> TF round-trip
[num_rt, den_rt] = ss2tf(A, B, C, D);
num_rt = num_rt / den_rt(1);
den_rt = den_rt / den_rt(1);

fprintf('pole0_re,pole0_im,pole1_re,pole1_im,D_00,num_rt_0,num_rt_1,den_rt_0,den_rt_1,den_rt_2\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        real(p(1)), imag(p(1)), real(p(2)), imag(p(2)), ...
        D(1,1), ...
        num_rt(1), num_rt(2), ...
        den_rt(1), den_rt(2), den_rt(3));
