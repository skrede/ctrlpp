# Trajectory Generation

Trajectory generation computes time-parameterised motion profiles that
guide a system from one configuration to another while respecting kinematic
and dynamic constraints. It is a fundamental component of motion control for
industrial machines, robots, and autonomous vehicles
[1, Ch. 1, pp. 1--14].

This page covers point-to-point motion planning, velocity profiles, polynomial
trajectories, spline interpolation, and online trajectory generation. The
primary reference is Biagiotti and Melchiorri (2009) [1], which provides a
comprehensive treatment of trajectory planning for automatic machines.

## Point-to-Point Motion

The simplest trajectory problem is moving from position $q_0$ to position
$q_1$ in time $T$, subject to constraints on velocity, acceleration, and
possibly jerk. The trajectory $q(t)$ must satisfy boundary conditions
[1, Ch. 3, pp. 49--72]:

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
at both endpoints [1, Sec. 3.2, pp. 53--56]:

$$
\begin{bmatrix} q_0 \\ q_1 \\ v_0 \\ v_1 \end{bmatrix} =
\begin{bmatrix} 1 & 0 & 0 & 0 \\ 1 & T & T^2 & T^3 \\ 0 & 1 & 0 & 0 \\ 0 & 1 & 2T & 3T^2 \end{bmatrix}
\begin{bmatrix} a_0 \\ a_1 \\ a_2 \\ a_3 \end{bmatrix}
$$

### Quintic Polynomials

A quintic polynomial additionally constrains accelerations at the endpoints,
providing smoother motion with continuous jerk
[1, Sec. 3.3, pp. 57--62]:

$$
q(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3 + a_4 t^4 + a_5 t^5
$$

Six coefficients satisfy six boundary conditions: position, velocity, and
acceleration at both $t = 0$ and $t = T$.

### Minimum-Jerk and Minimum-Snap

Higher-order polynomials minimise specific motion derivatives
[1, Sec. 3.4, pp. 63--66]:

- **Minimum-jerk** (5th order): minimises $\int_0^T \dddot{q}^2 \, dt$,
  producing smooth, human-like motions
- **Minimum-snap** (7th order): minimises $\int_0^T q^{(4)2} \, dt$,
  common in quadrotor trajectory planning where snap is proportional to
  motor force rate

## Velocity Profiles

Velocity profiles specify how the system accelerates, cruises, and
decelerates. They provide more intuitive control over the motion than
polynomial coefficients.

### Trapezoidal Velocity Profile

The trapezoidal profile has three phases [1, Sec. 3.5, pp. 67--80]:

1. **Acceleration**: constant acceleration $a_{\max}$ from rest to
   cruise velocity $v_{\max}$
2. **Cruise**: constant velocity $v_{\max}$
3. **Deceleration**: constant deceleration $-a_{\max}$ to rest

The position trajectory is piecewise quadratic-linear-quadratic. The time
durations of each phase are determined by the displacement $h = q_1 - q_0$,
maximum velocity $v_{\max}$, and maximum acceleration $a_{\max}$
[1, Sec. 3.5.1, pp. 68--72]:

$$
T_a = \frac{v_{\max}}{a_{\max}}, \quad T_v = \frac{h}{v_{\max}} - \frac{v_{\max}}{a_{\max}}, \quad T = T_v + 2 T_a
$$

If $h < v_{\max}^2 / a_{\max}$, the cruise phase vanishes and the profile
becomes triangular.

### Double-S (S-curve) Velocity Profile

The double-S profile limits jerk as well as acceleration, producing smoother
motion with reduced mechanical vibration. It has seven phases
[1, Sec. 3.6, pp. 80--106]:

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
and the displacement $h$ [1, Sec. 3.6.2, pp. 85--95].

Degenerate cases arise when the displacement is too small for all phases
to occur. The cruise phase may vanish, the constant-acceleration phases
may vanish, or both -- requiring careful case analysis
[1, Sec. 3.6.3, pp. 96--100].

## Spline Interpolation

When the trajectory must pass through multiple waypoints
$q_0, q_1, \ldots, q_n$ at times $t_0, t_1, \ldots, t_n$, spline
interpolation constructs a piecewise polynomial that is globally smooth
[1, Ch. 4, pp. 107--156].

### Cubic Splines

Cubic splines use third-order polynomials on each interval $[t_i, t_{i+1}]$,
with continuity of position, velocity, and acceleration at the knots. The
$n$ segments require $4n$ coefficients, determined by
[1, Sec. 4.2, pp. 112--120]:

- $2n$ interpolation conditions ($q(t_i) = q_i$ at each end of each segment)
- $n - 1$ velocity continuity conditions
- $n - 1$ acceleration continuity conditions
- $2$ boundary conditions (natural, clamped, or periodic)

This yields a tridiagonal linear system that is solved in $O(n)$ time.

### B-Splines

B-splines provide a basis function representation that supports local control:
moving one control point affects only a few nearby segments. The B-spline
curve of degree $p$ is [1, Sec. 4.4, pp. 130--140]:

$$
q(t) = \sum_{i=0}^{n} N_{i,p}(t) \, P_i
$$

where $N_{i,p}(t)$ are the B-spline basis functions (computed by the
Cox-de Boor recursion) and $P_i$ are the control points.

## Online Trajectory Generation

Online (real-time) trajectory generation computes trajectories on-the-fly in
response to changing targets, enabling reactive motion
[1, Ch. 7, pp. 237--280].

### Third-Order Online Planner

A third-order online planner limits velocity, acceleration, and jerk. At each
sample time, it evaluates the current kinematic state
$(q_k, \dot{q}_k, \ddot{q}_k)$ and computes the next sample toward the
target position, respecting all constraints
[1, Sec. 7.3, pp. 255--270]:

The planner constructs a sequence of constant-jerk phases:
1. Brake to zero acceleration if needed (to prepare for direction change)
2. Accelerate toward the target
3. Cruise at maximum velocity
4. Decelerate to arrive with zero velocity at the target

When the target changes mid-motion, the planner re-evaluates from the current
state, using a brake-to-zero-then-replan strategy for robustness.

## Multi-Axis Synchronisation

When multiple axes must reach their targets simultaneously (e.g., a multi-DOF
robot), the fastest axis determines the total motion time $T$, and slower
axes are scaled to match [1, Ch. 9, pp. 325--360]:

$$
T = \max_j T_j
$$

Each axis $j$ then adjusts its velocity and acceleration limits to complete
its motion in exactly $T$ while maintaining the same profile shape. The
scaling preserves the jerk, acceleration, and velocity constraints while
ensuring all axes start and stop together.

## References

[1] L. Biagiotti and C. Melchiorri, "Trajectory Planning for Automatic
Machines and Robots," Springer, 2009.

[2] N. S. Nise, "Control Systems Engineering," 7th ed., Wiley, 2015.

[3] K. M. Lynch and F. C. Park, "Modern Robotics: Mechanics, Planning, and
Control," Cambridge University Press, 2017.
