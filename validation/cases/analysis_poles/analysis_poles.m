% analysis_poles.m -- Reference system analysis via Octave Control Toolbox.
% Poles, stability, controllability, observability for continuous and discrete systems.
%
% Usage: octave --no-gui analysis_poles.m > analysis_poles_octave.csv

pkg load control;

% Continuous mass-spring-damper
A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = [1, 0];

sys_c = ss(A_c, B_c, C_c, 0);
p_c = pole(sys_c);

% Discrete version
sys_d = c2d(sys_c, 0.05);
p_d = pole(sys_d);

% Controllability rank
ctrb_rank = rank(ctrb(sys_d.A, sys_d.B));
obsv_rank = rank(obsv(sys_d.A, sys_d.C));

% Stability checks: 1.0 for stable, 0.0 for unstable
is_stable_c = all(real(p_c) < 0);
is_stable_d = all(abs(p_d) < 1);

fprintf('pc0_re,pc0_im,pc1_re,pc1_im,pd0_re,pd0_im,pd1_re,pd1_im,ctrb_rank,obsv_rank,stable_c,stable_d\n');

% Sort poles by imaginary part (positive first)
[~, idx_c] = sort(imag(p_c), 'descend');
p_c = p_c(idx_c);
[~, idx_d] = sort(imag(p_d), 'descend');
p_d = p_d(idx_d);

fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        real(p_c(1)), imag(p_c(1)), real(p_c(2)), imag(p_c(2)), ...
        real(p_d(1)), imag(p_d(1)), real(p_d(2)), imag(p_d(2)), ...
        ctrb_rank, obsv_rank, is_stable_c, is_stable_d);
