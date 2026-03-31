% lqr_step_response.m -- LQR closed-loop step response via Octave.
% Mass-spring-damper with LQR feedback, initial condition x0 = [1; 0].
%
% Usage: octave --no-gui lqr_step_response.m > lqr_step_response_octave.csv

pkg load control;

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = eye(2);
D_c = zeros(2, 1);

dt = 0.05;
duration = 10.0;

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, dt);

Ad = sys_d.A;
Bd = sys_d.B;

Q = [10 0; 0 1];
R = [1];

K = dlqr(Ad, Bd, Q, R);

x = [1; 0];

fprintf('time,x0,x1,u\n');

t = 0.0;
while t < duration - dt/2
    u = -K * x;
    fprintf('%.15e,%.15e,%.15e,%.15e\n', t, x(1), x(2), u);
    x = Ad * x + Bd * u;
    t = t + dt;
end
