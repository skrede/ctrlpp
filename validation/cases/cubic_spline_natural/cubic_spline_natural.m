% cubic_spline_natural.m -- Natural cubic spline via Octave built-in spline() + ppder().
%
% Usage: octave --no-gui cubic_spline_natural.m > cubic_spline_natural_octave.csv

% Waypoints
t_knots = [0.0, 1.0, 2.5, 4.0, 5.0];
q_knots = [0.0, 1.5, 0.8, 2.0, 1.0];

% Octave spline() uses not-a-knot BCs by default.
% For natural BCs (zero second derivative at endpoints), use csape from splines pkg.
pkg load splines;
pp = csape(t_knots, q_knots, 'variational');  % 'variational' = natural BCs

% Derivatives via ppder
pp_vel = ppder(pp);
pp_acc = ppder(pp_vel);

% Evaluate at fine grid
dt_eval = 0.05;

fprintf('time,position,velocity,acceleration\n');

t = t_knots(1);
while t <= t_knots(end) + dt_eval/2
    pos = ppval(pp, t);
    vel = ppval(pp_vel, t);
    acc = ppval(pp_acc, t);
    fprintf('%.15e,%.15e,%.15e,%.15e\n', t, pos, vel, acc);
    t = t + dt_eval;
end
