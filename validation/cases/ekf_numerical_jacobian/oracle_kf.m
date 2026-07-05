% oracle_kf.m -- Independent Octave oracle for the linear constant-velocity
% Kalman recursion exercised by the "ekf with numerical Jacobians converges
% on linear system" test case (tests/unit/ekf_convergence_test.cpp).
%
% Mirrors oracle_kf.py exactly (same recursion, Joseph-form covariance
% update matching ctrlpp/estimation/ekf.h::update_covariance) but
% implemented independently in Octave, so the two runtimes serve as
% mutually-independent oracles for the same closed-form linear KF math.
%
% Run: octave --no-gui oracle_kf.m

dt = 0.1;
F = [1.0, dt; 0.0, 1.0];
G = [0.5 * dt * dt; dt];
H = [1.0, 0.0];
Q = eye(2) * 0.01;
R = 1.0;

x = [0.0; 0.0];
P = eye(2) * 10.0;

true_pos = 0.0;
true_vel = 1.0;

for i = 0:49
    true_pos = true_pos + true_vel * dt;
    u = 0.0;

    % predict
    x = F * x + G * u;
    P = F * P * F' + Q;

    % measurement
    z = true_pos + 0.1 * sin(i);

    % update (Joseph form)
    z_pred = H * x;
    innovation = z - z_pred;
    S = H * P * H' + R;
    K = P * H' / S;
    x = x + K * innovation;
    IKH = eye(2) - K * H;
    P = IKH * P * IKH' + K * R * K';
end

printf("final_pos = %.17g\n", x(1));
printf("final_vel = %.17g\n", x(2));
printf("final_P00 = %.17g\n", P(1, 1));
printf("final_P11 = %.17g\n", P(2, 2));
printf("true_pos  = %.17g\n", true_pos);
printf("true_vel  = %.17g\n", true_vel);
