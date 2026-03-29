% pid_linear_step.m -- Linear PI controller output via Octave Control Toolbox.
% Constructs matching PI transfer function and uses lsim() for open-loop comparison.
%
% ctrlpp uses backward Euler integration: I(z) = Ki*Ts*z/(z-1)
% Matching TF: C(z) = Kp + Ki*Ts*z/(z-1) = ((Kp + Ki*Ts)*z - Kp) / (z - 1)
%
% Usage: octave --no-gui pid_linear_step.m > pid_linear_step_octave.csv

pkg load control;

Kp = 2.0;
Ki = 1.0;
Kd = 0.0;
dt = 0.01;

% Build PI TF matching ctrlpp backward Euler convention
Ki_dt = Ki * dt;
C = tf([Kp + Ki_dt, -Kp], [1, -1], dt);

% Known error signal: decaying sinusoid
n_steps = 500;
t_vec = (0:n_steps-1)' * dt;
e = exp(-0.5 * t_vec) .* sin(2 * pi * 0.5 * t_vec);

% Compute controller output via lsim
[u_out] = lsim(C, e, t_vec);

fprintf('time,error,control\n');
for k = 1:n_steps
    fprintf('%.15e,%.15e,%.15e\n', t_vec(k), e(k), u_out(k));
end
