% so3_exp_log.m -- SO(3) exp/log maps via Octave quaternion package.
%
% Usage: octave --no-gui so3_exp_log.m > so3_exp_log_octave.csv

pkg load quaternion;

% Test rotation vectors (angle-axis representation)
test_vecs = [
    0.0, 0.0, 0.0;
    0.1, 0.0, 0.0;
    0.0, 0.5, 0.0;
    0.0, 0.0, 1.0;
    1.0, 1.0, 1.0;
    pi*0.99, 0.0, 0.0;
    0.001, 0.002, 0.003;
];

n_tests = size(test_vecs, 1);

fprintf('phi_x,phi_y,phi_z,q_w,q_x,q_y,q_z,phi_rt_x,phi_rt_y,phi_rt_z\n');

for i = 1:n_tests
    phi = test_vecs(i, :)';
    theta = norm(phi);

    % Rotation vector -> quaternion via rot2q (axis-angle to quaternion)
    if theta < 1e-10
        q = quaternion(1, 0, 0, 0);
    else
        axis = phi / theta;
        q = rot2q(axis', theta);
    end

    % Extract components (w, x, y, z)
    q_w = q.w;
    q_x = q.x;
    q_y = q.y;
    q_z = q.z;

    % Canonicalize: ensure w >= 0
    if q_w < 0
        q_w = -q_w;
        q_x = -q_x;
        q_y = -q_y;
        q_z = -q_z;
    end

    % Quaternion -> rotation vector via q2rot
    q_canon = quaternion(q_w, q_x, q_y, q_z);
    [axis_rt, angle_rt] = q2rot(q_canon);
    phi_rt = angle_rt * axis_rt';

    % Handle identity case
    if abs(angle_rt) < 1e-10
        phi_rt = [0; 0; 0];
    end

    fprintf('%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n', ...
            phi(1), phi(2), phi(3), q_w, q_x, q_y, q_z, ...
            phi_rt(1), phi_rt(2), phi_rt(3));
end
