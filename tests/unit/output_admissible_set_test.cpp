#include "ctrlpp/mpc/terminal_set.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/mpc/qp_formulation.h"
#include "ctrlpp/mpc/output_admissible_set.h"

#include "ctrlpp/control/dare.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <variant>
#include <optional>
#include <algorithm>
#include <type_traits>

namespace
{

constexpr double dt = 0.1;
constexpr double inf = std::numeric_limits<double>::infinity();
// The QP encodes a polytope's absent lower side with this placeholder rather than -inf.
constexpr double absent_lower_side = -1e20;

template <std::size_t N>
using polytope = ctrlpp::polytopic_set<double, N>;

template <std::size_t N>
using row_vector = Eigen::Matrix<double, 1, static_cast<int>(N)>;

template <std::size_t N>
struct vertex
{
    ctrlpp::Vector<double, N> point;
    ctrlpp::Vector<double, N> error;
};

template <std::size_t N, std::size_t M>
struct plant
{
    ctrlpp::Matrix<double, N, N> A;
    ctrlpp::Matrix<double, N, M> B;
};

template <std::size_t N, std::size_t M>
struct box_bounds
{
    ctrlpp::Vector<double, M> u_min;
    ctrlpp::Vector<double, M> u_max;
    ctrlpp::Vector<double, N> x_min;
    ctrlpp::Vector<double, N> x_max;
};

plant<2, 1> double_integrator()
{
    return {.A = (ctrlpp::Matrix<double, 2, 2>() << 1.0, dt, 0.0, 1.0).finished(),
            .B = (ctrlpp::Matrix<double, 2, 1>() << 0.5 * dt * dt, dt).finished()};
}

template <std::size_t N, std::size_t M>
box_bounds<N, M> symmetric_bounds(double u_bound, double x_bound)
{
    return {.u_min = ctrlpp::Vector<double, M>::Constant(-u_bound),
            .u_max = ctrlpp::Vector<double, M>::Constant(u_bound),
            .x_min = ctrlpp::Vector<double, N>::Constant(-x_bound),
            .x_max = ctrlpp::Vector<double, N>::Constant(x_bound)};
}

template <std::size_t N, std::size_t M>
ctrlpp::Matrix<double, M, N> lqr_gain(const plant<N, M>& sys, const ctrlpp::Matrix<double, N, N>& P, const ctrlpp::Matrix<double, M, M>& R)
{
    ctrlpp::Matrix<double, M, M> const RpBtPB = R + sys.B.transpose() * P * sys.B;
    return -(RpBtPB.ldlt().solve(sys.B.transpose() * P * sys.A));
}

template <std::size_t N>
void append_row(polytope<N>& set, const row_vector<N>& row, double bound)
{
    auto const m = set.H.rows();
    set.H.conservativeResize(m + 1, Eigen::NoChange);
    set.h.conservativeResize(m + 1);
    set.H.row(m) = row;
    set.h(m) = bound;
}

// The residual G v - g is an (n + 1)-term inner product, so its computed value is within gamma_{n+1} of
// |G| |v| + |g| (Higham 2002, Sec. 3.1).
template <std::size_t N>
double residual_bound(const row_vector<N>& G, double g, const ctrlpp::Vector<double, N>& v)
{
    return ctrlpp::detail::rounding_bound<double>(static_cast<std::int32_t>(N + 1), G.cwiseAbs().dot(v.cwiseAbs()) + std::abs(g), 1.0);
}

template <std::size_t N>
double implication_bound(const row_vector<N>& G, double g, const vertex<N>& v)
{
    return residual_bound<N>(G, g, v.point) + G.cwiseAbs().dot(v.error);
}

// The vertex's forward error is bounded from the residual of its defining rows, |dv| <= |G_S^{-1}| (|r| + gamma_{n+1}
// (|G_S| |v| + |g_S|)) (Higham 2002, Sec. 7.2).
template <std::size_t N>
std::optional<vertex<N>> vertex_at(const polytope<N>& set, const std::vector<Eigen::Index>& subset)
{
    ctrlpp::Matrix<double, N, N> G;
    ctrlpp::Vector<double, N> g;
    for(std::size_t k = 0; k < N; ++k)
    {
        G.row(static_cast<Eigen::Index>(k)) = set.H.row(subset[k]);
        g(static_cast<Eigen::Index>(k)) = set.h(subset[k]);
    }
    Eigen::FullPivLU<ctrlpp::Matrix<double, N, N>> const lu(G);
    if(!lu.isInvertible())
        return std::nullopt;
    ctrlpp::Vector<double, N> const p = lu.solve(g);
    ctrlpp::Vector<double, N> slack;
    for(Eigen::Index k = 0; k < static_cast<Eigen::Index>(N); ++k)
        slack(k) = std::abs(G.row(k).dot(p) - g(k)) + residual_bound<N>(G.row(k), g(k), p);
    return vertex<N>{.point = p, .error = lu.inverse().cwiseAbs() * slack};
}

template <std::size_t N>
bool admitted(const polytope<N>& set, const std::vector<bool>& chosen, const vertex<N>& v)
{
    for(Eigen::Index j = 0; j < set.H.rows(); ++j)
        if(!chosen[static_cast<std::size_t>(j)] && set.H.row(j).dot(v.point) - set.h(j) > implication_bound<N>(set.H.row(j), set.h(j), v))
            return false;
    return true;
}

template <std::size_t N>
std::vector<vertex<N>> vertices_of(const polytope<N>& set)
{
    std::vector<vertex<N>> found;
    auto const m = static_cast<std::size_t>(set.H.rows());
    if(m < N)
        return found;
    std::vector<bool> chosen(m, false);
    std::fill_n(chosen.begin(), N, true);
    std::vector<Eigen::Index> subset;
    do
    {
        subset.clear();
        for(std::size_t i = 0; i < m; ++i)
            if(chosen[i])
                subset.push_back(static_cast<Eigen::Index>(i));
        if(auto v = vertex_at<N>(set, subset); v && admitted<N>(set, chosen, *v))
            found.push_back(*v);
    } while(std::prev_permutation(chosen.begin(), chosen.end()));
    return found;
}

// The rows the QP is posed with for the terminal stage, read back from the assembled constraint matrix and
// bounds: one row per finite upper side and one per present lower side.
template <std::size_t N, std::size_t M>
polytope<N> qp_terminal_rows(const ctrlpp::terminal_set<double, N>& tset, const plant<N, M>& sys, const box_bounds<N, M>& b)
{
    constexpr int horizon = 10;
    std::optional<ctrlpp::Vector<double, N>> const x_min{b.x_min}, x_max{b.x_max};
    std::optional<ctrlpp::Vector<double, M>> const u_min{b.u_min}, u_max{b.u_max}, du_max{};
    std::optional<ctrlpp::terminal_set<double, N>> const terminal{tset};
    auto const layout = ctrlpp::detail::constraint_layout<double, N, M>(horizon, x_min, x_max, u_min, u_max, du_max, terminal);
    Eigen::MatrixXd const matrix(ctrlpp::detail::build_constraint_matrix<double, N, M>(horizon, sys.A, sys.B, x_min, x_max, false, u_min, u_max, du_max, terminal));
    auto const [lower, upper] = ctrlpp::detail::build_bounds_vectors<double, N, M>(horizon, ctrlpp::Vector<double, N>::Zero(), x_min, x_max, u_min, u_max, du_max, terminal);
    polytope<N> rows{.H = Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>(0, N), .h = Eigen::VectorXd(0)};
    for(int r = layout.terminal_row; r < layout.n_con; ++r)
    {
        row_vector<N> const row = matrix.block(r, horizon * static_cast<int>(N), 1, static_cast<int>(N));
        if(std::isfinite(upper(r)))
            append_row<N>(rows, row, upper(r));
        if(lower(r) != absent_lower_side)
            append_row<N>(rows, -row, -lower(r));
    }
    return rows;
}

struct invariance_record
{
    double worst_ratio;
    double worst_margin;
    double worst_use;
};

// H_i (A_K v) - h_i is a two-product chain followed by a subtraction, within gamma_{2n+1} of |H_i| |A_K| |v| + |h_i|
// (Higham 2002, Sec. 3.5); the vertex's forward error reaches it through |H_i| |A_K| (Sec. 7.2). The construction may
// also have dropped the row H_i A_K as implied within its own residual bound, so that allowance is counted once more.
template <std::size_t N>
double invariance_bound(const row_vector<N>& H_i, double h_i, const ctrlpp::Matrix<double, N, N>& A_K, const vertex<N>& v)
{
    row_vector<N> const chain = H_i.cwiseAbs() * A_K.cwiseAbs();
    double const scale = chain.dot(v.point.cwiseAbs()) + std::abs(h_i);
    double const propagated = chain.dot(v.error);
    double const check = ctrlpp::detail::rounding_bound<double>(static_cast<std::int32_t>(2 * N + 1), scale, 1.0) + propagated;
    double const construction = ctrlpp::detail::rounding_bound<double>(static_cast<std::int32_t>(N + 1), scale, 1.0) + propagated;
    return check + construction;
}

template <std::size_t N>
invariance_record invariance_of(const polytope<N>& set, const ctrlpp::Matrix<double, N, N>& A_K)
{
    auto const vertices = vertices_of<N>(set);
    REQUIRE(vertices.size() >= N + 1);
    invariance_record record{.worst_ratio = -inf, .worst_margin = -inf, .worst_use = -inf};
    for(auto const& v : vertices)
    {
        ctrlpp::Vector<double, N> const image = A_K * v.point;
        for(Eigen::Index i = 0; i < set.H.rows(); ++i)
        {
            double const value = set.H.row(i).dot(image);
            double const bound = invariance_bound<N>(set.H.row(i), set.h(i), A_K, v);
            record.worst_ratio = std::max(record.worst_ratio, value / set.h(i));
            record.worst_margin = std::max(record.worst_margin, value - set.h(i) - bound);
            record.worst_use = std::max(record.worst_use, (value - set.h(i)) / bound);
        }
    }
    return record;
}

template <std::size_t N>
double worst_violation(const polytope<N>& rows, const std::vector<vertex<N>>& vertices)
{
    double worst = -inf;
    for(auto const& v : vertices)
        for(Eigen::Index i = 0; i < rows.H.rows(); ++i)
            worst = std::max(worst, rows.H.row(i).dot(v.point) - rows.h(i) - implication_bound<N>(rows.H.row(i), rows.h(i), v));
    return worst;
}

template <std::size_t N>
polytope<N> without_row(const polytope<N>& set, Eigen::Index dropped)
{
    polytope<N> reduced{.H = Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>(0, N), .h = Eigen::VectorXd(0)};
    for(Eigen::Index j = 0; j < set.H.rows(); ++j)
        if(j != dropped)
            append_row<N>(reduced, set.H.row(j), set.h(j));
    return reduced;
}

template <std::size_t N, std::size_t M>
polytope<N> output_rows(const ctrlpp::Matrix<double, M, N>& K, const box_bounds<N, M>& b)
{
    polytope<N> rows{.H = Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>(0, N), .h = Eigen::VectorXd(0)};
    auto add = [&rows](const row_vector<N>& row, double lower, double upper)
    {
        if(std::isfinite(upper))
            append_row<N>(rows, row, upper);
        if(std::isfinite(lower))
            append_row<N>(rows, -row, -lower);
    };
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
        add(row_vector<N>::Unit(j), b.x_min(j), b.x_max(j));
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(M); ++i)
        add(K.row(i), b.u_min(i), b.u_max(i));
    return rows;
}

// Every output row C A_K^k for k = 0..last, unit-normalized and unpruned.
template <std::size_t N>
polytope<N> horizon_rows(const polytope<N>& outputs, const ctrlpp::Matrix<double, N, N>& A_K, std::int64_t last)
{
    polytope<N> rows{.H = Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>(0, N), .h = Eigen::VectorXd(0)};
    ctrlpp::Matrix<double, N, N> power = ctrlpp::Matrix<double, N, N>::Identity();
    for(std::int64_t k = 0; k <= last; ++k, power = power * A_K)
        for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
        {
            row_vector<N> const row = outputs.H.row(i) * power;
            if(row.norm() > 0.0)
                append_row<N>(rows, row / row.norm(), outputs.h(i) / row.norm());
        }
    return rows;
}

// The Stein solution X = sum_k (A_K^T)^k A_K^k, summed until a term no longer changes X in floating point; an
// independent route to the matrix the construction solves for in Kronecker form.
template <std::size_t N>
ctrlpp::Matrix<double, N, N> lyapunov_series(const ctrlpp::Matrix<double, N, N>& A_K)
{
    ctrlpp::Matrix<double, N, N> X = ctrlpp::Matrix<double, N, N>::Identity();
    ctrlpp::Matrix<double, N, N> power = A_K;
    for(;;)
    {
        ctrlpp::Matrix<double, N, N> const next = X + power.transpose() * power;
        if(next == X)
            return X;
        X = next;
        power = power * A_K;
    }
}

// The a-priori step bound: V(x) = x^T X x contracts by 1 - 1 / lambda_max(X), r2 bounds V on O_{n-1} and r1 is the
// largest level of V inside the output constraints.
template <std::size_t N>
std::int64_t step_cap(const ctrlpp::Matrix<double, N, N>& A_K, const polytope<N>& outputs)
{
    auto const X = lyapunov_series<N>(A_K);
    ctrlpp::Matrix<double, N, N> const X_inv = X.inverse();
    double const lambda_max = Eigen::SelfAdjointEigenSolver<ctrlpp::Matrix<double, N, N>>(X).eigenvalues().maxCoeff();
    double r1 = inf;
    for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
        r1 = std::min(r1, outputs.h(i) * outputs.h(i) / (outputs.H.row(i) * X_inv * outputs.H.row(i).transpose())(0));
    double r2 = 0.0;
    for(auto const& v : vertices_of<N>(horizon_rows<N>(outputs, A_K, static_cast<std::int64_t>(N) - 1)))
        r2 = std::max(r2, v.point.dot(X * v.point));
    double const k_star = r2 <= r1 ? 0.0 : std::ceil(std::log(r1 / r2) / std::log1p(-1.0 / lambda_max));
    return std::max(static_cast<std::int64_t>(N), static_cast<std::int64_t>(k_star));
}

template <std::size_t N, std::size_t M>
struct lqr_case
{
    ctrlpp::terminal_ingredients_result<double, N, M> ingredients;
    ctrlpp::Matrix<double, M, N> K;
    ctrlpp::Matrix<double, N, N> A_K;
};

template <std::size_t N, std::size_t M>
lqr_case<N, M> solve_case(const plant<N, M>& sys, const ctrlpp::Matrix<double, M, M>& R, const box_bounds<N, M>& b)
{
    auto const ingredients = ctrlpp::terminal_ingredients<double, N, M>(sys.A, sys.B, ctrlpp::Matrix<double, N, N>::Identity(), R, b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(ingredients.has_value());
    auto const K = lqr_gain<N, M>(sys, ingredients->Qf, R);
    return {.ingredients = *ingredients, .K = K, .A_K = sys.A + sys.B * K};
}

ctrlpp::Matrix<double, 1, 1> scalar(double value)
{
    return ctrlpp::Matrix<double, 1, 1>::Constant(value);
}

template <typename Refusal, typename Result>
std::optional<Refusal> refusal_as(const Result& result)
{
    if(result.has_value())
        return std::nullopt;
    if constexpr(std::is_same_v<typename Result::error_type, Refusal>)
        return result.error();
    else if(auto const* cause = std::get_if<Refusal>(&result.error()))
        return *cause;
    return std::nullopt;
}

// The pose is proved by the Riccati solve on the same operands, and terminal_ingredients must carry that solve's own cause.
void require_riccati_cause(const plant<2, 1>& sys, const ctrlpp::Matrix<double, 2, 2>& Q, ctrlpp::dare_error cause)
{
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    CHECK(refusal_as<ctrlpp::dare_error>(ctrlpp::dare<double, 2, 1>(sys.A, sys.B, Q, scalar(1.0))) == cause);
    auto const ingredients = ctrlpp::terminal_ingredients<double, 2, 1>(sys.A, sys.B, Q, scalar(1.0), b.u_min, b.u_max, b.x_min, b.x_max);
    CHECK_FALSE(ingredients.has_value());
    CHECK(refusal_as<ctrlpp::dare_error>(ingredients) == cause);
}

plant<2, 1> diagonal_plant(double a1, double a2, double b1, double b2)
{
    return {.A = ctrlpp::Vector<double, 2>(a1, a2).asDiagonal(), .B = (ctrlpp::Matrix<double, 2, 1>() << b1, b2).finished()};
}

}

TEST_CASE("The ellipsoid option's QP box is not invariant under the LQR law", "[terminal_set][output_admissible]")
{
    auto const sys = double_integrator();
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    auto const solved = solve_case<2, 1>(sys, scalar(1.0), b);

    auto const box = qp_terminal_rows<2, 1>(ctrlpp::terminal_set<double, 2>{solved.ingredients.set}, sys, b);
    REQUIRE(box.H.rows() == 4);

    auto const record = invariance_of<2>(box, solved.A_K);
    INFO(std::setprecision(17) << "worst one-step ratio " << record.worst_ratio << ", margin over the bound " << record.worst_margin);
    CHECK(record.worst_margin > 0.0);
}

TEST_CASE("The output admissible set the QP receives is invariant under the LQR law", "[terminal_set][output_admissible]")
{
    auto const sys = double_integrator();
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    auto const solved = solve_case<2, 1>(sys, scalar(1.0), b);
    REQUIRE(solved.ingredients.invariant_set.has_value());
    auto const& set = *solved.ingredients.invariant_set;

    auto const posed = qp_terminal_rows<2, 1>(ctrlpp::terminal_set<double, 2>{set}, sys, b);
    // One posed row per QP row means every lower side carried the absent-side placeholder.
    REQUIRE(posed.H.rows() == set.H.rows());
    CHECK(posed.H == set.H);
    CHECK(posed.h == set.h);
    // The minimal representation of a full-dimensional polytope is its facet list, which an independent
    // linear-programming construction of this set counts at 30.
    CHECK(set.H.rows() == 30);

    auto const record = invariance_of<2>(posed, solved.A_K);
    INFO(std::setprecision(17) << "rows " << set.H.rows() << ", worst ratio " << record.worst_ratio << ", largest use of the bound " << record.worst_use);
    CHECK(record.worst_margin <= 0.0);
}

TEST_CASE("Every row of the output admissible set is a facet", "[terminal_set][output_admissible]")
{
    auto const solved = solve_case<2, 1>(double_integrator(), scalar(1.0), symmetric_bounds<2, 1>(1.0, 10.0));
    REQUIRE(solved.ingredients.invariant_set.has_value());
    auto const& set = *solved.ingredients.invariant_set;

    for(Eigen::Index i = 0; i < set.H.rows(); ++i)
    {
        INFO("row " << i);
        polytope<2> const row{.H = set.H.row(i), .h = set.h.segment(i, 1)};
        CHECK(worst_violation<2>(row, vertices_of<2>(without_row<2>(set, i))) > 0.0);
    }
}

TEST_CASE("Without state bounds the output admissible set is bounded and invariant", "[terminal_set][output_admissible]")
{
    auto const sys = double_integrator();
    auto const b = symmetric_bounds<2, 1>(1.0, inf);
    auto const solved = solve_case<2, 1>(sys, scalar(0.1), b);
    REQUIRE(solved.ingredients.invariant_set.has_value());
    auto const& set = *solved.ingredients.invariant_set;

    auto const vertices = vertices_of<2>(set);
    REQUIRE(vertices.size() >= 3);
    for(auto const& v : vertices)
        CHECK(v.point.allFinite());
    // An independent vertex-enumeration construction of this set counts 16 facets.
    CHECK(set.H.rows() == 16);

    auto const posed = qp_terminal_rows<2, 1>(ctrlpp::terminal_set<double, 2>{set}, sys, b);
    REQUIRE(posed.H.rows() == set.H.rows());
    auto const record = invariance_of<2>(posed, solved.A_K);
    INFO(std::setprecision(17) << "rows " << set.H.rows() << ", worst ratio " << record.worst_ratio << ", largest use of the bound " << record.worst_use);
    CHECK(record.worst_margin <= 0.0);
}

TEST_CASE("The output admissible set agrees with every output row up to the a-priori step bound", "[terminal_set][output_admissible]")
{
    auto const x_bound = GENERATE(10.0, inf);
    auto const R = x_bound == inf ? scalar(0.1) : scalar(1.0);
    auto const b = symmetric_bounds<2, 1>(1.0, x_bound);
    auto const solved = solve_case<2, 1>(double_integrator(), R, b);
    REQUIRE(solved.ingredients.invariant_set.has_value());
    auto const& set = *solved.ingredients.invariant_set;

    auto const outputs = output_rows<2, 1>(solved.K, b);
    auto const cap = step_cap<2>(solved.A_K, outputs);
    auto const horizon = horizon_rows<2>(outputs, solved.A_K, cap);
    INFO("state bound " << x_bound << ", step bound " << cap << ", horizon rows " << horizon.H.rows());
    CHECK(worst_violation<2>(horizon, vertices_of<2>(set)) <= 0.0);
    CHECK(worst_violation<2>(set, vertices_of<2>(horizon)) <= 0.0);
}

TEST_CASE("The output admissible set refuses each failed hypothesis by name", "[terminal_set][output_admissible]")
{
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    ctrlpp::Matrix<double, 1, 2> const K = (ctrlpp::Matrix<double, 1, 2>() << -1.0, -1.5).finished();
    ctrlpp::Matrix<double, 2, 2> const stable = (ctrlpp::Matrix<double, 2, 2>() << 0.9, 0.1, 0.0, 0.8).finished();

    SECTION("an unstable closed loop")
    {
        ctrlpp::Matrix<double, 2, 2> const unstable = 1.1 * ctrlpp::Matrix<double, 2, 2>::Identity();
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(unstable, K, b.u_min, b.u_max, b.x_min, b.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::closed_loop_not_stable);
    }
    SECTION("a non-finite gain")
    {
        ctrlpp::Matrix<double, 1, 2> const broken = (ctrlpp::Matrix<double, 1, 2>() << -1.0, std::numeric_limits<double>::quiet_NaN()).finished();
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(stable, broken, b.u_min, b.u_max, b.x_min, b.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::non_finite_input);
    }
    SECTION("an infinite closed-loop entry")
    {
        ctrlpp::Matrix<double, 2, 2> broken = stable;
        broken(0, 1) = inf;
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(broken, K, b.u_min, b.u_max, b.x_min, b.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::non_finite_input);
    }
    SECTION("a NaN bound with an unstable closed loop is reported as non-finite data first")
    {
        ctrlpp::Matrix<double, 2, 2> const unstable = 1.1 * ctrlpp::Matrix<double, 2, 2>::Identity();
        ctrlpp::Vector<double, 2> const x_min(std::numeric_limits<double>::quiet_NaN(), -10.0);
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(unstable, K, b.u_min, b.u_max, x_min, b.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::non_finite_input);
    }
    SECTION("an infinite bound side is absent, not corrupt")
    {
        ctrlpp::Vector<double, 2> const x_min(-inf, -10.0);
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(stable, K, b.u_min, b.u_max, x_min, b.x_max);
        CHECK(result.has_value());
    }
    SECTION("a finite state bound that does not straddle zero")
    {
        ctrlpp::Vector<double, 2> const x_min(0.5, -10.0);
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(stable, K, b.u_min, b.u_max, x_min, b.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::origin_not_interior);
    }
    SECTION("a NaN bound")
    {
        ctrlpp::Vector<double, 2> const x_max(10.0, std::numeric_limits<double>::quiet_NaN());
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(stable, K, b.u_min, b.u_max, b.x_min, x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::non_finite_input);
    }
    SECTION("finitely bounded outputs that do not observe the closed loop")
    {
        ctrlpp::Matrix<double, 2, 2> const diagonal = ctrlpp::Vector<double, 2>(0.5, 0.8).asDiagonal();
        ctrlpp::Matrix<double, 1, 2> const first_only = (ctrlpp::Matrix<double, 1, 2>() << 1.0, 0.0).finished();
        auto const unbounded = symmetric_bounds<2, 1>(1.0, inf);
        auto const result = ctrlpp::compute_output_admissible_set<double, 2, 1>(diagonal, first_only, unbounded.u_min, unbounded.u_max, unbounded.x_min, unbounded.x_max);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == ctrlpp::terminal_set_error::unbounded_terminal_set);
    }
}

TEST_CASE("The step budget refuses only below the last step the set needs, and its default changes nothing", "[terminal_set][output_admissible]")
{
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    auto const solved = solve_case<2, 1>(double_integrator(), scalar(1.0), b);
    auto const within = [&](std::int64_t budget)
    { return ctrlpp::compute_output_admissible_set<double, 2, 1>(solved.A_K, solved.K, b.u_min, b.u_max, b.x_min, b.x_max, 1.0, budget); };

    auto const by_default = ctrlpp::compute_output_admissible_set<double, 2, 1>(solved.A_K, solved.K, b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(by_default.has_value());
    CHECK(by_default->H.rows() == 30);
    auto const unbounded = within(std::numeric_limits<std::int64_t>::max());
    REQUIRE(unbounded.has_value());
    CHECK(unbounded->H == by_default->H);
    CHECK(unbounded->h == by_default->h);

    auto const cap = step_cap<2>(solved.A_K, output_rows<2, 1>(solved.K, b));
    std::int64_t smallest = -1;
    for(std::int64_t budget = 0; budget <= cap && smallest < 0; ++budget)
    {
        auto const result = within(budget);
        if(result.has_value())
            smallest = budget;
        else
            CHECK(result.error() == ctrlpp::terminal_set_error::not_converged);
    }
    INFO("smallest admitted budget " << smallest << ", step bound " << cap);
    REQUIRE(smallest > 0);
    auto const at = within(smallest);
    REQUIRE(at.has_value());
    CHECK(at->H == by_default->H);
    CHECK(at->h == by_default->h);
    auto const below = within(smallest - 1);
    REQUIRE_FALSE(below.has_value());
    CHECK(below.error() == ctrlpp::terminal_set_error::not_converged);
}

TEST_CASE("terminal_ingredients reports an unsupported dimension through the invariant set alone", "[terminal_set][output_admissible]")
{
    ctrlpp::Matrix<double, 5, 5> const A = 0.5 * ctrlpp::Matrix<double, 5, 5>::Identity();
    ctrlpp::Matrix<double, 5, 1> const B = ctrlpp::Matrix<double, 5, 1>::Ones();
    ctrlpp::Matrix<double, 5, 5> const Q = ctrlpp::Matrix<double, 5, 5>::Identity();
    auto const b = symmetric_bounds<5, 1>(1.0, inf);

    auto const ingredients = ctrlpp::terminal_ingredients<double, 5, 1>(A, B, Q, scalar(1.0), b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(ingredients.has_value());
    REQUIRE_FALSE(ingredients->invariant_set.has_value());
    CHECK(ingredients->invariant_set.error() == ctrlpp::terminal_set_error::dimension_unsupported);

    auto const P = ctrlpp::dare<double, 5, 1>(A, B, Q, scalar(1.0));
    REQUIRE(P.has_value());
    CHECK(ingredients->Qf == P->P);
    auto const K = lqr_gain<5, 1>(plant<5, 1>{.A = A, .B = B}, P->P, scalar(1.0));
    auto const ellipsoid = ctrlpp::compute_ellipsoidal_set<double, 5, 1>(P->P, K, b.u_min, b.u_max, b.x_min, b.x_max);
    REQUIRE(ellipsoid.has_value());
    CHECK(ingredients->set.P == ellipsoid->P);
    CHECK(ingredients->set.alpha == ellipsoid->alpha);
}

TEST_CASE("terminal_ingredients refuses a singular input weight whatever A is", "[terminal_set][singular_r]")
{
    auto const b = symmetric_bounds<2, 1>(1.0, 10.0);
    ctrlpp::Matrix<double, 2, 2> const Q = ctrlpp::Matrix<double, 2, 2>::Identity();
    plant<2, 1> nilpotent{.A = (ctrlpp::Matrix<double, 2, 2>() << 0.0, 1.0, 0.0, 0.0).finished(), .B = (ctrlpp::Matrix<double, 2, 1>() << 0.0, 1.0).finished()};
    auto ingredients_of = [&](const plant<2, 1>& sys, double r) { return ctrlpp::terminal_ingredients<double, 2, 1>(sys.A, sys.B, Q, scalar(r), b.u_min, b.u_max, b.x_min, b.x_max); };

    SECTION("a zero input weight")
    {
        CHECK(refusal_as<ctrlpp::terminal_set_error>(ingredients_of(double_integrator(), 0.0)) == ctrlpp::terminal_set_error::singular_input_weight);
    }

    SECTION("a rank-deficient A with a definite input weight is the Riccati solve's own refusal")
    {
        CHECK(refusal_as<ctrlpp::dare_error>(ingredients_of(nilpotent, 1.0)) == ctrlpp::dare_error::singular_a);
    }

    SECTION("a rank-deficient A does not hide a zero input weight")
    {
        CHECK(refusal_as<ctrlpp::terminal_set_error>(ingredients_of(nilpotent, 0.0)) == ctrlpp::terminal_set_error::singular_input_weight);
    }
}

TEST_CASE("terminal_ingredients carries the Riccati solve's own refusal", "[terminal_set][riccati_refusal]")
{
    ctrlpp::Matrix<double, 2, 2> const Q = ctrlpp::Matrix<double, 2, 2>::Identity();

    SECTION("singular_a: a nilpotent A")
    {
        require_riccati_cause({.A = (ctrlpp::Matrix<double, 2, 2>() << 0.0, 1.0, 0.0, 0.0).finished(), .B = (ctrlpp::Matrix<double, 2, 1>() << 0.0, 1.0).finished()}, Q,
                              ctrlpp::dare_error::singular_a);
    }

    SECTION("non_stabilizable: a zero state weight on a Jordan block at one")
    {
        require_riccati_cause(double_integrator(), ctrlpp::Matrix<double, 2, 2>::Zero(), ctrlpp::dare_error::non_stabilizable);
    }

    SECTION("arithmetic_limit: a controllable pair whose solution keeps less than half the significand")
    {
        require_riccati_cause(diagonal_plant(2.0, 0.5, 1e-6, 1.0), Q, ctrlpp::dare_error::arithmetic_limit);
    }

    SECTION("singular_u11: an exactly uncontrollable unstable mode")
    {
        require_riccati_cause(diagonal_plant(2.0, 0.5, 0.0, 1.0), Q, ctrlpp::dare_error::singular_u11);
    }

    SECTION("an ellipsoid refusal stays a terminal-set refusal")
    {
        auto b = symmetric_bounds<2, 1>(1.0, 10.0);
        b.u_min = scalar(0.5);
        auto const sys = double_integrator();
        auto const ingredients = ctrlpp::terminal_ingredients<double, 2, 1>(sys.A, sys.B, Q, scalar(1.0), b.u_min, b.u_max, b.x_min, b.x_max);
        CHECK(refusal_as<ctrlpp::terminal_set_error>(ingredients) == ctrlpp::terminal_set_error::input_zero_not_interior);
    }
}

TEST_CASE("The output admissible set is constructed and invariant at four states", "[terminal_set][output_admissible]")
{
    plant<4, 2> sys{.A = ctrlpp::Matrix<double, 4, 4>::Zero(), .B = ctrlpp::Matrix<double, 4, 2>::Zero()};
    auto const single = double_integrator();
    for(Eigen::Index axis = 0; axis < 2; ++axis)
    {
        sys.A.block<2, 2>(2 * axis, 2 * axis) = single.A;
        sys.B.block<2, 1>(2 * axis, axis) = single.B;
    }
    auto const b = symmetric_bounds<4, 2>(1.0, 10.0);
    auto const solved = solve_case<4, 2>(sys, ctrlpp::Matrix<double, 2, 2>::Identity(), b);
    REQUIRE(solved.ingredients.invariant_set.has_value());
    auto const& set = *solved.ingredients.invariant_set;

    auto const posed = qp_terminal_rows<4, 2>(ctrlpp::terminal_set<double, 4>{set}, sys, b);
    REQUIRE(posed.H.rows() == set.H.rows());
    auto const record = invariance_of<4>(posed, solved.A_K);
    INFO(std::setprecision(17) << "rows " << set.H.rows() << ", worst ratio " << record.worst_ratio << ", largest use of the bound " << record.worst_use);
    CHECK(record.worst_margin <= 0.0);
}
