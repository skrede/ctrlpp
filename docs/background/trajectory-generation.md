# Trajectory Generation

Trajectory generation computes time-parametrized motion profiles that
guide a system from one configuration to another while respecting kinematic
and dynamic constraints. It is a fundamental component of motion control for
industrial machines and robots [1, Ch. 1, pp. 1--12].

This page covers point-to-point motion planning, velocity profiles, polynomial
trajectories, spline interpolation, and online trajectory generation. The
primary reference is Biagiotti and Melchiorri (2008) [1], which provides a
comprehensive treatment of trajectory planning for automatic machines.

## Point-to-Point Motion

The simplest trajectory problem is moving from position $q_0$ to position
$q_1$ in time $T$, subject to constraints on velocity, acceleration, and
possibly jerk. The trajectory $q(t)$ must satisfy boundary conditions
[1, Sec. 2.1, p. 15]:

$$
q(0) = q_0, \quad q(T) = q_1
$$

$$
\dot{q}(0) = v_0, \quad \dot{q}(T) = v_1
$$

Additional constraints may include:

- Maximum velocity: $|\dot{q}(t)| \le v_{\max}$
- Maximum acceleration: $|\ddot{q}(t)| \le a_{\max}$
- Maximum jerk: $|\dddot{q}(t)| \le j_{\max}$

## Polynomial Trajectories

### Cubic Polynomials

A cubic polynomial $q(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3$ has four
coefficients, sufficient to satisfy position and velocity boundary conditions
at both endpoints [1, Sec. 2.1.4, eq. (2.1)--(2.2), pp. 23--24]:

$$
\begin{bmatrix} q_0 \\ q_1 \\ v_0 \\ v_1 \end{bmatrix} =
\begin{bmatrix} 1 & 0 & 0 & 0 \\ 1 & T & T^2 & T^3 \\ 0 & 1 & 0 & 0 \\ 0 & 1 & 2T & 3T^2 \end{bmatrix}
\begin{bmatrix} a_0 \\ a_1 \\ a_2 \\ a_3 \end{bmatrix}
$$

ctrlpp writes the polynomial in normalized time $\tau = t / T$, with
coefficients $c_k = a_k T^k$, and evaluates it by Horner's method. With
$h = q_1 - q_0$ the solution of the system above is

$$
c_0 = q_0, \quad c_1 = v_0 T, \quad c_2 = 3h - T (2 v_0 + v_1), \quad c_3 = -2h + T (v_0 + v_1)
$$

### Quintic Polynomials

A quintic polynomial additionally constrains accelerations at the endpoints,
so a segment joins its neighbors with continuous acceleration when the
boundary accelerations match theirs
[1, Sec. 2.1.5, eq. (2.4)--(2.5), pp. 26--27]:

$$
q(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + a_4 t^4 + a_5 t^5
$$

Six coefficients satisfy six boundary conditions: position, velocity, and
acceleration at both $t = 0$ and $t = T$. A septic polynomial adds the jerk
at both ends, for eight conditions and eight coefficients.

### Minimum-Jerk and Minimum-Snap

Higher-order polynomials minimize specific motion derivatives:

- **Minimum-jerk** (5th order): among all motions meeting the position,
  velocity and acceleration conditions at both ends, the quintic minimizes
  $\int_0^T \dddot{q}^2 \, dt$
- **Minimum-snap** (7th order): with the jerk conditions added, the septic
  minimizes $\int_0^T (q^{(4)})^2 \, dt$

ctrlpp builds these polynomials from the boundary conditions; it runs no
optimization.

## Velocity Profiles

Velocity profiles specify how the system accelerates, cruises, and
decelerates. They provide more intuitive control over the motion than
polynomial coefficients.

### Trapezoidal Velocity Profile

For a rest-to-rest move with $h > 0$ the trapezoidal profile has three
phases [1, Sec. 3.2, pp. 62--64]:

1. **Acceleration**: constant acceleration $a_{\max}$ from rest to
   cruise velocity $v_{\max}$
2. **Cruise**: constant velocity $v_{\max}$
3. **Deceleration**: constant deceleration $-a_{\max}$ to rest

The position trajectory is piecewise quadratic-linear-quadratic. The time
durations of each phase are determined by the displacement $h = q_1 - q_0$,
maximum velocity $v_{\max}$, and maximum acceleration $a_{\max}$
[1, Sec. 3.2.2, eq. (3.8), pp. 65--66]:

$$
T_a = \frac{v_{\max}}{a_{\max}}, \quad T_v = \frac{h}{v_{\max}} - \frac{v_{\max}}{a_{\max}}, \quad T = T_v + 2 T_a
$$

If $h < v_{\max}^2 / a_{\max}$, the cruise phase vanishes and the profile
becomes triangular. ctrlpp's trapezoidal profile also takes nonzero boundary
velocities $v_0$ and $v_1$, for which the phase durations follow from the
same three-phase structure.

### Double-S (S-curve) Velocity Profile

The double-S profile limits jerk as well as acceleration, so the acceleration
is continuous. It has seven phases
[1, Sec. 3.4, pp. 79--80]:

1. Increasing acceleration (jerk $= j_{\max}$)
2. Constant acceleration ($a_{\max}$)
3. Decreasing acceleration (jerk $= -j_{\max}$)
4. Cruise (constant velocity)
5. Increasing deceleration (jerk $= -j_{\max}$)
6. Constant deceleration ($-a_{\max}$)
7. Decreasing deceleration (jerk $= j_{\max}$)

The position trajectory is piecewise: cubic during jerk phases, quadratic
during constant-acceleration phases, and linear during cruise. The phase
durations depend on the constraints $v_{\max}$, $a_{\max}$, $j_{\max}$
and the displacement $h$ [1, Sec. 3.4, eq. (3.17)--(3.27), pp. 80--84].

Degenerate cases arise when the displacement is too small for all phases
to occur. The cruise phase may vanish, the constant-acceleration phases
may vanish, or both &mdash; requiring careful case analysis
[1, Sec. 3.4, pp. 81--84; Sec. 3.4.3, pp. 90--93].

## Spline Interpolation

When the trajectory must pass through multiple waypoints
$q_0, q_1, \ldots, q_n$ at times $t_0, t_1, \ldots, t_n$, spline
interpolation constructs a piecewise polynomial that is globally smooth
[1, Ch. 4, pp. 151--219].

### Cubic Splines

Cubic splines use third-order polynomials on each interval $[t_i, t_{i+1}]$,
with continuity of position, velocity, and acceleration at the knots. The
$n$ segments require $4n$ coefficients, determined by
[1, Sec. 4.4, pp. 166--167]:

- $2n$ interpolation conditions ($q(t_i) = q_i$ at each end of each segment)
- $n - 1$ velocity continuity conditions
- $n - 1$ acceleration continuity conditions
- $2$ boundary conditions (natural, clamped, or periodic)

ctrlpp takes the velocities at the knots as the unknowns. Natural and
clamped end conditions give a tridiagonal linear system, solved in $O(n)$
time by the Thomas algorithm; periodic end conditions give a cyclic
tridiagonal one, also solved in $O(n)$.

### Smoothing Splines

A smoothing spline trades interpolation for smoothness when the waypoints
are noisy. ctrlpp's smoothing spline has zero acceleration at both ends and
solves for the second derivatives $d$ at the interior knots:

$$
(R + \lambda Q^\top Q) \, d = Q^\top q, \qquad \lambda = \frac{2 (1 - \mu)}{3 \mu}
$$

where $R$ is the tridiagonal smoothness matrix of the knot spacings, $Q$
the second-difference matrix and $q$ the waypoints. The smoothed knot
positions are $q - \lambda Q d$. The weight $\mu \in (0, 1]$ sets the
tradeoff: $\mu = 1$ gives $\lambda = 0$ and the interpolating spline, and
as $\mu$ approaches 0 the curve approaches the least-squares straight line.

### B-Splines

B-splines provide a basis function representation that supports local control:
moving one control point affects only a few nearby segments. The B-spline
curve of degree $p$ is [1, Sec. 4.5, eq. (4.40)--(4.41), pp. 194--195]:

$$
q(t) = \sum_{i=0}^{n} N_{i,p}(t) \, P_i
$$

where $N_{i,p}(t)$ are the B-spline basis functions (computed by the
Cox-de Boor recursion) and $P_i$ are the control points, scalars in ctrlpp.
ctrlpp evaluates the curve, and its derivatives from derivative control
points, by de Boor's algorithm rather than by summing basis functions.

## Online Trajectory Generation

Online (real-time) trajectory generation computes trajectories on-the-fly in
response to changing targets, enabling reactive motion
[1, Sec. 4.6, pp. 208--209].

### Third-Order Online Planner

A third-order online planner limits velocity, acceleration, and jerk
[1, Sec. 4.6.1, p. 209]. When the target changes, it plans a double-S
profile from the current kinematic state $(q, \dot{q}, \ddot{q})$ to the
target at rest; between target changes each sample evaluates that plan.

The plan is a sequence of constant-jerk phases:
1. Bring a nonzero starting acceleration to zero
2. Accelerate toward the target
3. Cruise at maximum velocity
4. Decelerate to arrive with zero velocity at the target

A move in the direction the axis is already traveling carries the current
velocity into the profile. A velocity pointing away from the target, or too
large to stop in the remaining distance, is braked to rest first, and the
profile is planned from the stopping point. The second-order planner does
the same with velocity and acceleration limits only.

## Multi-Axis Synchronization

When multiple axes must reach their targets simultaneously (e.g., a multi-DOF
robot), the slowest axis determines the total motion time $T$, and the
other axes are slowed to match [1, Sec. 3.2.3, pp. 66--67; Sec. 5.3,
pp. 241--243]:

$$
T = \max_j T_j
$$

Each axis $j$ is then retimed to complete its motion in $T$. ctrlpp's
trapezoidal profile lowers its cruise velocity at a fixed acceleration; the
double-S profile is rebuilt under velocity, acceleration and jerk limits
scaled by the first, second and third power of one factor. Every axis is
checked before any is retimed, so a set that cannot be synchronized is left
unchanged.

## References

[1] L. Biagiotti and C. Melchiorri, "Trajectory Planning for Automatic
Machines and Robots," Springer, 2008.

[2] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.

[3] K. M. Lynch and F. C. Park, "Modern Robotics: Mechanics, Planning, and
Control," Cambridge University Press, 2017.
