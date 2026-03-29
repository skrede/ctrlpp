% lqi_tracking.m -- LQI (LQR with integral action) step tracking.
% Mass-spring-damper with position output, tracks step reference.
%
% Usage: octave --no-gui lqi_tracking.m > lqi_tracking_octave.csv

pkg load control;

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = [1, 0];
D_c = [0];
dt = 0.05;
duration = 15.0;

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, dt);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

nx = 2;
ny = 1;
nu = 1;
nx_aug = nx + ny;

% Augmented system for integral action
A_aug = zeros(nx_aug);
A_aug(1:nx, 1:nx) = Ad;
A_aug(nx+1:end, 1:nx) = -Cd;
A_aug(nx+1:end, nx+1:end) = eye(ny);

B_aug = zeros(nx_aug, nu);
B_aug(1:nx, :) = Bd;

% Weight matrices for augmented system
Q_aug = diag([10, 1, 50]);  % position, velocity, integral
R_aug = [0.1];

K_aug = dlqr(A_aug, B_aug, Q_aug, R_aug);
Kx = K_aug(1:nu, 1:nx);
Ki = K_aug(1:nu, nx+1:end);

x = [0; 0];
xi = 0;  % integral state
ref = 1.0;

fprintf('time,position,velocity,integral,control\n');

t = 0.0;
while t < duration - dt/2
    y = Cd * x;
    u = -Kx * x - Ki * xi;

    fprintf('%.15e,%.15e,%.15e,%.15e,%.15e\n', t, x(1), x(2), xi, u);

    x = Ad * x + Bd * u;
    xi = xi + (ref - y);

    t = t + dt;
end
