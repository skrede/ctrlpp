#ifndef HPP_GUARD_CTRLPP_TESTS_ESTIMATION_NOMALLOC_GRID_H
#define HPP_GUARD_CTRLPP_TESTS_ESTIMATION_NOMALLOC_GRID_H

// Shared machinery for the estimator rows' no-allocation dimension grids: the
// corpus the published stack figures were taken on, one builder per row, the
// warm-up-then-arm helper, and the compile-time point list a row's case walks.
//
// nomalloc_harness.h MUST be included before this header, because this one
// pulls in Eigen. That ordering is not left to convention: the harness errors
// out if Eigen was parsed before it.
//
// WHY THE GRID LIVES IN SEVERAL TRANSLATION UNITS. Every point is a distinct
// instantiation chain, and the compile cost of a point is dominated by that
// instantiation rather than by the dimension it instantiates at. One
// translation unit holding every row's grid is therefore a single-core
// serialization point that no build parallelism can reduce, and its peak
// resident set is the sum of every row's instantiations in one process. One
// unit per row is the same total work spread over as many cores as the build
// has, at the cost of one repeated header parse per unit, one more link, and
// one more serialized test. The harness's own contract already forces one
// executable per unit, so the split costs no additional serialization beyond
// that.
//
// WHAT EACH ROW ARMS. Every configuration the real-time safety matrix prints a
// per-configuration figure for, which is three groups:
//
//   1. the row's measurement table -- five rungs on each of three lines through
//      its grid, with the rung where the lines cross counted once;
//   2. every value its supported-maximum table names as fitting a task stack,
//      which is what a caller reads off that document and instantiates;
//   3. the largest value each axis instantiates at, and the corner where both
//      axes are largest;
//
// plus both sides of the boundary below, which the grid found and which is
// published beside the figures it bounds.
//
// THE HEAP CLAIM HAS A BOUNDARY, AND IT IS AT A MEASUREMENT DIMENSION OF 48.
// The gain solve on four of these rows is a column-pivoting Householder QR of
// the NY x NY innovation covariance applied to a multi-column right-hand side.
// Eigen applies a Householder sequence BY BLOCK once the sequence is at least
// BlockSize = 48 long and the destination has more than one column
// (HouseholderSequence.h), and the blocked application declares its block
// reflector factor as
//
//     Matrix<Scalar, TFactorSize, TFactorSize, RowMajor> T(nbVecs, nbVecs)
//
// with TFactorSize taken from the destination's compile-time column count
// (BlockHouseholder.h). The destination there is a dynamically sized block, so
// TFactorSize is Dynamic and T is HEAP ALLOCATED -- one 48 x 48 block of 18,432
// bytes per update, through Eigen's own aligned allocator.
//
// So the 48 below is read out of the linear-algebra library rather than fitted,
// and it is walked from both sides here: 47 is allocation-free and 48 is not.
// Two exemptions are measured rather than assumed -- a single-state filter,
// whose right-hand side is one column, stays on the unblocked path at any
// measurement dimension; and this row's ukf, whose gain decomposition defaults
// to LDLT, has no Householder sequence at all until a caller selects the QR
// option, which reintroduces the boundary exactly.
//
// THE GLOBAL COUNTER CANNOT SEE THAT ALLOCATION. It goes through Eigen's
// aligned_malloc, which calls std::malloc directly, so the counting operator
// new reads zero at every allocating point below. A guard carrying only the
// counter would report every one of them allocation-free, which is why the
// allocating points assert BOTH that the sentinel fires and that the counter
// stayed at zero.
//
// WHAT IS NOT ARMED, stated because a zero across a grid reads as a claim about
// everything between its points: the interior fill behind the published
// supported maxima is 4,677 points and is not individually armed. Nothing here
// is a claim about a dimension lying between two armed ones.

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/manifold_ukf.h"
#include "ctrlpp/estimation/estimation_types.h"

#include "ctrlpp/lie/so3.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <utility>


namespace ctrlpp_test::grid
{

using ctrlpp::Vector;
using ctrlpp::Matrix;

// Runs the callable inside an armed no-malloc window and returns the number of
// heap allocations it performed. The count is sampled before the guard is
// released and before any test macro runs, so framework-internal allocations
// cannot pollute it.
template <typename Fn>
std::size_t guarded_allocations(Fn&& fn)
{
    ctrlpp_test::scoped_no_malloc guard;
    std::forward<Fn>(fn)();
    return guard.allocations();
}

// Steady-state predict/update pairs run inside a grid point's armed window.
// This is not a resolution parameter: the guard traps on the FIRST offending
// call, and every call after the warm-up runs the same code over the same
// storage. The repetition covers only a branch that would be taken on alternate
// steps, and it is kept small because the grid reaches dimensions where one
// pair is a 128 x 128 factorization. The single fixed-instantiation case each
// row also carries keeps its own 128-step window.
constexpr int steady_state_steps = 4;

// The input dimension every vector-state grid point is driven at, which is the
// one the published stack figures were taken at.
constexpr std::size_t input_dimension = 1;

// Warm-up-then-arm, shared by every grid point: one predict/update pair OUTSIDE
// the window flushes lazy one-time instantiation, then the same pair runs
// inside the armed window. The count is sampled inside guarded_allocations,
// before any Catch2 macro runs in the window, because framework macros may
// themselves allocate. BOTH mechanisms are asserted -- zero counted allocations
// and a clean Eigen-side sentinel -- because either alone is a silent false
// pass.
template <typename Filter, typename Input, typename Measurement>
void require_alloc_free_steady_state(Filter& filter, const Input& u, const Measurement& z)
{
    filter.predict(u);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    const std::size_t allocations = guarded_allocations([&] {
        for(int step = 0; step < steady_state_steps; ++step)
        {
            filter.predict(u);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE_FALSE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

// The other side of the boundary, armed the same way and asserted exactly: at
// and above a measurement dimension of 48 with a multi-column right-hand side,
// the gain solve allocates its block reflector factor on the heap once per
// update. The counter is required to stay at ZERO here, because that is the
// measured demonstration that it cannot see this allocation -- a guard carrying
// only the counter would report these configurations allocation-free.
template <typename Filter, typename Input, typename Measurement>
void require_library_allocation_steady_state(Filter& filter, const Input& u, const Measurement& z)
{
    filter.predict(u);
    REQUIRE(filter.update(z).has_value());

    bool all_stepped = true;
    const std::size_t allocations = guarded_allocations([&] {
        for(int step = 0; step < steady_state_steps; ++step)
        {
            filter.predict(u);
            all_stepped = all_stepped && filter.update(z).has_value();
        }
    });
    REQUIRE(ctrlpp_test::eigen_violation());
    REQUIRE(all_stepped);

    REQUIRE(allocations == 0);
    REQUIRE(filter.state().allFinite());
}

// One armed grid point, carrying the two dimensions of the row it belongs to.
// A row with a single caller axis puts its STRUCTURAL dimension in the second
// slot rather than sweeping it, so the point still reads as a configuration
// rather than as half of one.
template <std::size_t First, std::size_t Second>
struct point
{
    static constexpr std::size_t first  = First;
    static constexpr std::size_t second = Second;
};

// A row's grid. arm() applies the callable to every point and returns how many
// points it armed, so the count a case asserts on is produced by the traversal
// rather than restated beside it.
template <typename... Points>
struct over
{
    static constexpr std::size_t size = sizeof...(Points);

    template <typename Arm>
    static std::size_t arm(Arm&& arm_point)
    {
        (arm_point(Points{}), ...);
        return size;
    }
};


// The corpus, in the shape the stack-watermark instrument builds it, so an
// armed point and a published stack figure are the same configuration of the
// same plant rather than two different systems at equal dimensions. The three
// vector-state rows take the forward-Euler damped chain; the two attitude rows
// propagate a constant body rate on SO(3) and measure the gravity direction
// repeated to the output dimension. The two families are not comparable at
// equal state dimension because they are not the same plant.
template <std::size_t NX>
auto state_matrix() -> Matrix<double, NX, NX>
{
    constexpr double dt = 0.01;

    Matrix<double, NX, NX> A = Matrix<double, NX, NX>::Identity();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt * 1.0;
    return A;
}

// One input per group of states, acting on the last state of its group. The
// group size is clamped at one so the map stays defined when there are more
// inputs than states; at input dimension 1 the clamp is inactive.
template <std::size_t NX, std::size_t NU>
auto input_matrix() -> Matrix<double, NX, NU>
{
    constexpr double dt = 0.01;

    Matrix<double, NX, NU> B = Matrix<double, NX, NU>::Zero();
    const std::size_t      group = NU <= NX ? NX / NU : std::size_t{1};
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t reach = (j + 1) * group;
        B(int((reach < NX ? reach : NX) - 1), int(j)) = dt;
    }
    return B;
}

// Output i reads state i mod NX. The modulus is what makes the map defined at
// every pair of dimensions, including the points where the output dimension
// runs past the state dimension; there it duplicates rows, which leaves the
// innovation covariance nonsingular because the measurement noise is identity.
template <std::size_t NX, std::size_t NY>
auto output_matrix() -> Matrix<double, NY, NX>
{
    Matrix<double, NY, NX> H = Matrix<double, NY, NX>::Zero();
    for(std::size_t i = 0; i < NY; ++i)
        H(int(i), int(i % NX)) = 1.0;
    return H;
}

template <std::size_t NX, std::size_t NU>
struct chain_dynamics
{
    Matrix<double, NX, NX> F;
    Matrix<double, NX, NU> G;

    auto operator()(const Vector<double, NX>& x, const Vector<double, NU>& u) const -> Vector<double, NX>
    {
        return Vector<double, NX>{F * x + G * u};
    }

    auto jacobian_x(const Vector<double, NX>&, const Vector<double, NU>&) const -> Matrix<double, NX, NX>
    {
        return F;
    }

    auto jacobian_u(const Vector<double, NX>&, const Vector<double, NU>&) const -> Matrix<double, NX, NU>
    {
        return G;
    }
};

template <std::size_t NX, std::size_t NY>
struct chain_measurement
{
    Matrix<double, NY, NX> H;

    auto operator()(const Vector<double, NX>& x) const -> Vector<double, NY>
    {
        return Vector<double, NY>{H * x};
    }

    auto jacobian(const Vector<double, NX>&) const -> Matrix<double, NY, NX>
    {
        return H;
    }
};

struct body_rate_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, 3>& omega) const
        -> Eigen::Quaternion<double>
    {
        const Vector<double, 3> increment = (omega * dt).eval();
        return (q * ctrlpp::so3::exp(increment)).normalized();
    }
};

// The gravity direction in body axes, repeated to whatever output dimension is
// asked for. Repetition rather than truncation, so the map is defined above
// three outputs as well as below it.
template <std::size_t NY>
auto gravity_direction(const Eigen::Quaternion<double>& q) -> Vector<double, NY>
{
    const Matrix<double, 3, 3> rotation = q.toRotationMatrix();
    const Vector<double, 3>    down     = rotation.transpose().col(2);

    Vector<double, NY> z;
    for(std::size_t i = 0; i < NY; ++i)
        z(int(i)) = down(int(i % 3));
    return z;
}

template <std::size_t NY>
struct attitude_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q) const -> Vector<double, NY>
    {
        return gravity_direction<NY>(q);
    }
};

template <std::size_t NB, std::size_t NY>
struct attitude_bias_measurement
{
    auto operator()(const Eigen::Quaternion<double>& q, const Vector<double, NB>&) const -> Vector<double, NY>
    {
        return gravity_direction<NY>(q);
    }
};

// The constant body rate every attitude grid point is driven at.
inline auto body_rate() -> Vector<double, 3>
{
    Vector<double, 3> rate;
    rate << 0.01, -0.02, 0.03;
    return rate;
}

// The zero input every vector-state grid point is driven with.
inline auto zero_input() -> Vector<double, input_dimension>
{
    return Vector<double, input_dimension>::Zero();
}

template <std::size_t NX, std::size_t NY>
auto build_kalman()
{
    ctrlpp::discrete_state_space<double, NX, input_dimension, NY> system;
    system.A = state_matrix<NX>();
    system.B = input_matrix<NX, input_dimension>();
    system.C = output_matrix<NX, NY>();
    system.D = Matrix<double, NY, input_dimension>::Zero();

    const ctrlpp::kalman_config<double, NX, input_dimension, NY> config{
        .Q  = Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = Matrix<double, NY, NY>::Identity(),
        .x0 = Vector<double, NX>::Zero(),
        .P0 = Matrix<double, NX, NX>::Identity() * 10.0};

    return ctrlpp::kalman_filter<double, NX, input_dimension, NY>::create(system, config);
}

template <std::size_t NX, std::size_t NY>
auto build_ekf()
{
    const chain_dynamics<NX, input_dimension> dynamics{state_matrix<NX>(), input_matrix<NX, input_dimension>()};
    const chain_measurement<NX, NY>           measurement_map{output_matrix<NX, NY>()};

    const ctrlpp::ekf_config<double, NX, input_dimension, NY> config{
        .Q  = Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = Matrix<double, NY, NY>::Identity(),
        .x0 = Vector<double, NX>::Zero(),
        .P0 = Matrix<double, NX, NX>::Identity() * 10.0};

    return ctrlpp::ekf<double, NX, input_dimension, NY, chain_dynamics<NX, input_dimension>,
                       chain_measurement<NX, NY>>::create(dynamics, measurement_map, config);
}

template <std::size_t NX, std::size_t NY>
auto build_ukf()
{
    const chain_dynamics<NX, input_dimension> dynamics{state_matrix<NX>(), input_matrix<NX, input_dimension>()};
    const chain_measurement<NX, NY>           measurement_map{output_matrix<NX, NY>()};

    const ctrlpp::ukf_config<double, NX, input_dimension, NY> config{
        .Q  = Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = Matrix<double, NY, NY>::Identity(),
        .x0 = Vector<double, NX>::Zero(),
        .P0 = Matrix<double, NX, NX>::Identity() * 10.0};

    return ctrlpp::ukf<double, NX, input_dimension, NY, chain_dynamics<NX, input_dimension>,
                       chain_measurement<NX, NY>>::create(dynamics, measurement_map, config);
}

// The same row with its non-default gain decomposition selected, which is a
// public option and a different code path.
template <std::size_t NX, std::size_t NY>
auto build_ukf_qr()
{
    const chain_dynamics<NX, input_dimension> dynamics{state_matrix<NX>(), input_matrix<NX, input_dimension>()};
    const chain_measurement<NX, NY>           measurement_map{output_matrix<NX, NY>()};

    ctrlpp::ukf_config<double, NX, input_dimension, NY> config{
        .Q  = Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = Matrix<double, NY, NY>::Identity(),
        .x0 = Vector<double, NX>::Zero(),
        .P0 = Matrix<double, NX, NX>::Identity() * 10.0};
    config.decomposition = ctrlpp::gain_decomposition::qr;

    return ctrlpp::ukf<double, NX, input_dimension, NY, chain_dynamics<NX, input_dimension>,
                       chain_measurement<NX, NY>>::create(dynamics, measurement_map, config);
}

template <std::size_t NY>
auto build_manifold_ukf()
{
    ctrlpp::manifold_ukf_config<double, NY> config;
    config.Q *= 1e-6;

    return ctrlpp::manifold_ukf<double, NY, body_rate_dynamics, attitude_measurement<NY>>::create(
        body_rate_dynamics{}, attitude_measurement<NY>{}, config);
}

template <std::size_t NB, std::size_t NY>
auto build_mekf()
{
    ctrlpp::mekf_config<double, NB, NY> config;
    config.Q *= 1e-6;

    return ctrlpp::mekf<double, NB, NY, attitude_bias_measurement<NB, NY>>::create(
        attitude_bias_measurement<NB, NY>{}, config);
}

}

#endif
