% pole_placement.m -- Reference pole placement via Octave Control Toolbox.
% Single-input mass-spring-damper, place closed-loop poles at desired locations.
%
% Usage: octave --no-gui pole_placement.m > pole_placement_octave.csv

pkg load control;

A = [0, 1; -1, -0.5];
B = [0; 1];
C = eye(2);
D = zeros(2, 1);

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, 0.05);

Ad = sys_d.A;
Bd = sys_d.B;

% Desired poles inside unit circle
desired = [0.5 + 0.1i, 0.5 - 0.1i];
K = place(Ad, Bd, desired);

% Verify closed-loop eigenvalues
cl_eig = eig(Ad - Bd * K);

fprintf('K_00,K_01,cl_eig0_re,cl_eig0_im,cl_eig1_re,cl_eig1_im\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        K(1,1), K(1,2), real(cl_eig(1)), imag(cl_eig(1)), real(cl_eig(2)), imag(cl_eig(2)));
