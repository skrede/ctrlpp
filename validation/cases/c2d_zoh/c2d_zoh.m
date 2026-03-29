% c2d_zoh.m -- Reference ZOH discretisation via Octave Control Toolbox.
% Outputs flattened A_d, B_d matrices.
%
% Usage: octave --no-gui c2d_zoh.m > c2d_zoh_octave.csv

pkg load control;

% Mass-spring-damper: m=1, k=1, b=0.5
A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = eye(2);
D_c = zeros(2, 1);

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, 0.05);

Ad = sys_d.A;
Bd = sys_d.B;

fprintf('Ad_00,Ad_01,Ad_10,Ad_11,Bd_00,Bd_10\n');
fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
        Ad(1,1), Ad(1,2), Ad(2,1), Ad(2,2), Bd(1,1), Bd(2,1));
