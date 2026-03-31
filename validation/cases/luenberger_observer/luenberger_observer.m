% luenberger_observer.m -- Luenberger observer via Octave Control Toolbox.
% Uses place() for observer gain, estim() to build observer, lsim() for simulation.
%
% Usage: octave --no-gui luenberger_observer.m > luenberger_observer_octave.csv

pkg load control;

A_c = [0, 1; -1, -0.5];
B_c = [0; 1];
C_c = [1, 0];
D_c = [0];
dt = 0.05;
duration = 5.0;

sys_c = ss(A_c, B_c, C_c, D_c);
sys_d = c2d(sys_c, dt);
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;

% Observer gain via place on dual system (toolbox function)
desired_obs = [0.3, 0.2];
L = place(Ad', Cd', desired_obs)';

% Simulate observer manually using toolbox-computed gain L
% (estim+lsim would work for linear sim, but ctrlpp uses predict/update split)
x_true = [1; 0];
x_est = [0; 0];

fprintf('time,x_true_0,x_true_1,x_est_0,x_est_1\n');

t = 0.0;
while t < duration - dt/2
    u = [0.5 * sin(t)];

    fprintf('%.15e,%.15e,%.15e,%.15e,%.15e\n', t, x_true(1), x_true(2), x_est(1), x_est(2));

    % Predict
    x_est = Ad * x_est + Bd * u;
    x_true = Ad * x_true + Bd * u;

    % Update with toolbox-computed gain
    z_new = Cd * x_true;
    x_est = x_est + L * (z_new - Cd * x_est);

    t = t + dt;
end
