#ifndef HPP_GUARD_CTRLPP_TESTS_OUTPUT_ADMISSIBLE_ORACLE_H
#define HPP_GUARD_CTRLPP_TESTS_OUTPUT_ADMISSIBLE_ORACLE_H

// An admissibility oracle for the maximal output admissible set that shares no code with the library's polytope
// helpers: boundedness by non-negative least squares, admissibility by simulation to a Lyapunov certificate, and
// minimality of the returned rows by sampling.

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/detail/covariance_ops.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <algorithm>

namespace ctrlpp_test::oracle
{

template <std::size_t N>
using point = ctrlpp::Vector<double, N>;

template <std::size_t N>
using square = ctrlpp::Matrix<double, N, N>;

template <std::size_t N>
using row_vector = Eigen::Matrix<double, 1, static_cast<int>(N)>;

template <std::size_t N>
using rows = ctrlpp::polytopic_set<double, N>;

// The standard fixes mt19937's sequence but not its distributions', so every draw is formed from the raw output.
inline double unit_draw(std::mt19937& g)
{
    return static_cast<double>(g()) / (static_cast<double>(std::mt19937::max()) + 1.0);
}

inline double half_integer_draw(std::mt19937& g)
{
    return static_cast<double>(static_cast<std::int32_t>(g() % 5U) - 2) / 2.0;
}

inline double bound(std::int32_t count, double scale)
{
    return ctrlpp::detail::rounding_bound<double>(count, scale, 1.0);
}

template <std::size_t N>
rows<N> empty_rows()
{
    return {.H = Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>(0, static_cast<Eigen::Index>(N)), .h = Eigen::VectorXd(0)};
}

template <std::size_t N>
void append(rows<N>& set, const row_vector<N>& row, double h)
{
    auto const m = set.H.rows();
    set.H.conservativeResize(m + 1, Eigen::NoChange);
    set.h.conservativeResize(m + 1);
    set.H.row(m) = row;
    set.h(m) = h;
}

// One row per finite bound side; a zero row constrains nothing, since the origin is interior.
template <std::size_t N, std::size_t M>
rows<N> output_rows(const ctrlpp::Matrix<double, M, N>& K, const point<M>& u_min, const point<M>& u_max, const point<N>& x_min, const point<N>& x_max)
{
    rows<N> outputs = empty_rows<N>();
    auto add = [&outputs](const row_vector<N>& row, double lower, double upper)
    {
        if(row.norm() == 0.0)
            return;
        if(std::isfinite(upper))
            append<N>(outputs, row, upper);
        if(std::isfinite(lower))
            append<N>(outputs, -row, -lower);
    };
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
        add(row_vector<N>::Unit(j), x_min(j), x_max(j));
    for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(M); ++i)
        add(K.row(i), u_min(i), u_max(i));
    return outputs;
}

// The output rows C A_K^k for k = 0..last, unit-normalized; a zero row constrains nothing and is left out.
template <std::size_t N>
rows<N> horizon_rows(const rows<N>& outputs, const square<N>& A_K, std::int64_t last)
{
    rows<N> result = empty_rows<N>();
    square<N> power = square<N>::Identity();
    for(std::int64_t k = 0; k <= last; ++k, power = power * A_K)
        for(Eigen::Index i = 0; i < outputs.H.rows(); ++i)
        {
            row_vector<N> const row = outputs.H.row(i) * power;
            if(row.norm() > 0.0)
                append<N>(result, row / row.norm(), outputs.h(i) / row.norm());
        }
    return result;
}

// V(x) = x^T X x with X the Stein solution sum_k (A_K^T)^k A_K^k, summed until a term no longer changes X. With
// D = X - A_K^T X A_K positive definite, V(A_K x) <= (1 - decay) V(x) for decay = lambda_min(D) / lambda_max(X).
template <std::size_t N>
struct lyapunov_certificate
{
    square<N> X;
    square<N> X_inv;
    double lambda_min;
    double lambda_max;
    double decay;
};

template <std::size_t N>
lyapunov_certificate<N> certify(const square<N>& A_K)
{
    square<N> X = square<N>::Identity();
    for(square<N> power = A_K;; power = power * A_K)
    {
        square<N> const next = X + power.transpose() * power;
        if(next == X)
            break;
        X = next;
    }
    square<N> const D = X - A_K.transpose() * X * A_K;
    auto const spectrum = Eigen::SelfAdjointEigenSolver<square<N>>(X, Eigen::EigenvaluesOnly).eigenvalues();
    double const d_min = Eigen::SelfAdjointEigenSolver<square<N>>(D, Eigen::EigenvaluesOnly).eigenvalues().minCoeff();
    return {.X = X, .X_inv = X.inverse(), .lambda_min = spectrum.minCoeff(), .lambda_max = spectrum.maxCoeff(), .decay = d_min / spectrum.maxCoeff()};
}

template <std::size_t N>
struct closed_loop
{
    square<N> A_K;
    rows<N> outputs;
    lyapunov_certificate<N> lyapunov;
};

template <std::size_t N, std::size_t M>
closed_loop<N> make_closed_loop(const square<N>& A_K, const ctrlpp::Matrix<double, M, N>& K, const point<M>& u_min, const point<M>& u_max, const point<N>& x_min,
                                const point<N>& x_max)
{
    return {.A_K = A_K, .outputs = output_rows<N, M>(K, u_min, u_max, x_min, x_max), .lyapunov = certify<N>(A_K)};
}

// The largest level of V inside the outputs shifted by `offset`: max of c x over {V <= r} is sqrt(r c X^{-1} c^T).
template <std::size_t N>
double level(const closed_loop<N>& loop, double offset)
{
    double r = std::numeric_limits<double>::infinity();
    for(Eigen::Index i = 0; i < loop.outputs.H.rows(); ++i)
    {
        double const shifted = loop.outputs.h(i) + offset;
        r = std::min(r, shifted * shifted / (loop.outputs.H.row(i) * loop.lyapunov.X_inv * loop.outputs.H.row(i).transpose())(0));
    }
    return r;
}

// The steps the contraction needs to take V from v0 down to r, plus one for the rounding of the count itself.
inline std::int64_t certificate_steps(double r, double v0, double decay)
{
    if(v0 <= r)
        return 0;
    return static_cast<std::int64_t>(std::ceil(std::log(r / v0) / std::log1p(-decay))) + 1;
}

// x stays in the outputs shifted by `offset` for every k: no row is exceeded before V(x_k) <= r(offset), and from
// there the trajectory stays in that level set, which lies inside the shifted outputs.
template <std::size_t N>
bool admissible_forever(point<N> x, const closed_loop<N>& loop, double offset)
{
    double const r = level<N>(loop, offset);
    auto const steps = certificate_steps(r, x.dot(loop.lyapunov.X * x), loop.lyapunov.decay);
    for(std::int64_t k = 0; k <= steps; ++k, x = loop.A_K * x)
    {
        if(((loop.outputs.H * x - loop.outputs.h).array() > offset).any())
            return false;
        if(x.dot(loop.lyapunov.X * x) <= r)
            return true;
    }
    return false;
}

template <std::size_t N>
using generators = Eigen::Matrix<double, static_cast<int>(N), Eigen::Dynamic>;

template <std::size_t N>
Eigen::VectorXd passive_least_squares(const generators<N>& E, const point<N>& f, const std::vector<bool>& passive)
{
    std::vector<Eigen::Index> chosen;
    for(Eigen::Index j = 0; j < E.cols(); ++j)
        if(passive[static_cast<std::size_t>(j)])
            chosen.push_back(j);
    Eigen::MatrixXd const E_P = E(Eigen::all, chosen);
    Eigen::VectorXd const z_P = E_P.colPivHouseholderQr().solve(f);
    Eigen::VectorXd z = Eigen::VectorXd::Zero(E.cols());
    for(std::size_t k = 0; k < chosen.size(); ++k)
        z(chosen[k]) = z_P(static_cast<Eigen::Index>(k));
    return z;
}

// The free, unrejected index whose gradient w_j exceeds the rounding of its own evaluation the most, or -1 at optimality.
template <std::size_t N>
Eigen::Index entering_index(const generators<N>& E, const point<N>& f, const Eigen::VectorXd& x, const std::vector<bool>& excluded)
{
    auto const m = static_cast<std::int32_t>(E.cols());
    point<N> const r = f - E * x;
    point<N> const r_error = (E.cwiseAbs() * x.cwiseAbs() + f.cwiseAbs()).unaryExpr([m](double s) { return bound(m + 1, s); });
    Eigen::Index best = -1;
    double best_excess = 0.0;
    for(Eigen::Index j = 0; j < E.cols(); ++j)
    {
        double const excess = E.col(j).dot(r) - bound(static_cast<std::int32_t>(N), E.col(j).cwiseAbs().dot(r.cwiseAbs())) - E.col(j).cwiseAbs().dot(r_error);
        if(!excluded[static_cast<std::size_t>(j)] && excess > best_excess)
        {
            best = j;
            best_excess = excess;
        }
    }
    return best;
}

// Moves x toward z until the first passive coefficient reaches zero, which leaves the passive set.
inline void step_to_feasible(Eigen::VectorXd& x, const Eigen::VectorXd& z, std::vector<bool>& passive)
{
    double alpha = 1.0;
    Eigen::Index blocking = 0;
    for(Eigen::Index j = 0; j < x.size(); ++j)
    {
        double const ratio = x(j) > z(j) ? x(j) / (x(j) - z(j)) : 0.0;
        if(passive[static_cast<std::size_t>(j)] && z(j) <= 0.0 && ratio < alpha)
        {
            alpha = ratio;
            blocking = j;
        }
    }
    x += alpha * (z - x);
    x(blocking) = 0.0;
    for(Eigen::Index j = 0; j < x.size(); ++j)
        if(x(j) <= 0.0)
        {
            x(j) = 0.0;
            passive[static_cast<std::size_t>(j)] = false;
        }
}

// Lawson and Hanson 1974, Ch. 23, Sec. 3, algorithm NNLS. As in their code, the outer iterations are capped at 3 m, and
// a candidate whose coefficient would not come out positive is rejected until x next changes.
template <std::size_t N>
double nnls_residual(const generators<N>& E, const point<N>& f)
{
    auto const m = static_cast<std::size_t>(E.cols());
    Eigen::VectorXd x = Eigen::VectorXd::Zero(E.cols());
    std::vector<bool> passive(m, false);
    std::vector<bool> excluded(m, false);
    for(std::size_t iteration = 0; iteration < 3 * m; ++iteration)
    {
        Eigen::Index const t = entering_index<N>(E, f, x, excluded);
        if(t < 0)
            break;
        passive[static_cast<std::size_t>(t)] = true;
        if(passive_least_squares<N>(E, f, passive)(t) <= 0.0)
        {
            passive[static_cast<std::size_t>(t)] = false;
            excluded[static_cast<std::size_t>(t)] = true;
            continue;
        }
        for(std::size_t inner = 0; inner <= m; ++inner)
        {
            Eigen::VectorXd const z = passive_least_squares<N>(E, f, passive);
            bool interior = true;
            for(Eigen::Index j = 0; j < E.cols(); ++j)
                interior = interior && (!passive[static_cast<std::size_t>(j)] || z(j) > 0.0);
            if(interior)
            {
                x = z;
                break;
            }
            step_to_feasible(x, z, passive);
        }
        excluded = passive;
    }
    return (f - E * x).norm();
}

// The rows of H positively span R^N exactly when every +e_j and -e_j is a non-negative combination of them. If they do
// not, the cone they generate lies in a half-space {y : d^T y <= 0} with ||d|| = 1, and the target sign(d_j) e_j with
// |d_j| >= 1 / sqrt(N) is at least that far from the cone. The exact residual is therefore 0 or at least 1 / sqrt(N),
// and the decision is taken halfway between the two.
template <std::size_t N>
bool positively_spans(const Eigen::Matrix<double, Eigen::Dynamic, static_cast<int>(N)>& H)
{
    generators<N> const E = H.transpose();
    double const halfway = 0.5 / std::sqrt(static_cast<double>(N));
    for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
        for(double const sign : {1.0, -1.0})
            if(nnls_residual<N>(E, sign * point<N>::Unit(j)) >= halfway)
                return false;
    return true;
}

template <std::size_t N>
struct vertex
{
    point<N> p;
    point<N> error;
};

// The residual G v - g is an (N + 1)-term inner product within bound(N + 1, |G| |v| + |g|) (Higham 2002, Sec. 3.1);
// the vertex's forward error is |G_S^{-1}| times the residual of its defining rows (Higham 2002, Sec. 7.2).
template <std::size_t N>
double row_slack(const row_vector<N>& row, double h, const vertex<N>& v)
{
    return bound(static_cast<std::int32_t>(N + 1), row.cwiseAbs().dot(v.p.cwiseAbs()) + std::abs(h)) + row.cwiseAbs().dot(v.error);
}

template <std::size_t N>
std::optional<vertex<N>> solve_subset(const rows<N>& set, const std::vector<Eigen::Index>& subset)
{
    square<N> const G = set.H(subset, Eigen::all);
    point<N> const g = set.h(subset);
    Eigen::FullPivLU<square<N>> const lu(G);
    if(!lu.isInvertible())
        return std::nullopt;
    vertex<N> v{.p = lu.solve(g), .error = point<N>::Zero()};
    point<N> slack;
    for(Eigen::Index k = 0; k < static_cast<Eigen::Index>(N); ++k)
        slack(k) = std::abs(G.row(k).dot(v.p) - g(k)) + bound(static_cast<std::int32_t>(N + 1), G.row(k).cwiseAbs().dot(v.p.cwiseAbs()) + std::abs(g(k)));
    v.error = lu.inverse().cwiseAbs() * slack;
    return v;
}

template <std::size_t N>
bool satisfies_all(const rows<N>& set, const vertex<N>& v)
{
    for(Eigen::Index j = 0; j < set.H.rows(); ++j)
        if(set.H.row(j).dot(v.p) - set.h(j) > row_slack<N>(set.H.row(j), set.h(j), v))
            return false;
    return true;
}

// Every N-row subset is solved and kept when its point satisfies every row within that row's slack. The list places
// samples and is checked for admissibility; it never decides whether a row is redundant.
template <std::size_t N>
std::vector<vertex<N>> vertices_of(const rows<N>& set)
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
        if(auto const v = solve_subset<N>(set, subset); v && satisfies_all<N>(set, *v))
            found.push_back(*v);
    } while(std::prev_permutation(chosen.begin(), chosen.end()));
    return found;
}

template <std::size_t N>
struct sample_box
{
    point<N> lower;
    point<N> upper;
};

template <std::size_t N>
sample_box<N> doubled_box(const std::vector<vertex<N>>& vertices)
{
    sample_box<N> box{.lower = point<N>::Zero(), .upper = point<N>::Zero()};
    for(auto const& v : vertices)
    {
        box.lower = box.lower.cwiseMin(2.0 * v.p);
        box.upper = box.upper.cwiseMax(2.0 * v.p);
    }
    return box;
}

// The largest V over the box's corners and the vertices of the rows of the first N steps, which contain the set.
template <std::size_t N>
double largest_level(const closed_loop<N>& loop, const sample_box<N>& box)
{
    double v_max = 0.0;
    for(std::uint32_t corner = 0; corner < (1U << N); ++corner)
    {
        point<N> x;
        for(std::size_t j = 0; j < N; ++j)
            x(static_cast<Eigen::Index>(j)) = (corner >> j) & 1U ? box.upper(static_cast<Eigen::Index>(j)) : box.lower(static_cast<Eigen::Index>(j));
        v_max = std::max(v_max, x.dot(loop.lyapunov.X * x));
    }
    for(auto const& v : vertices_of<N>(horizon_rows<N>(loop.outputs, loop.A_K, static_cast<std::int64_t>(N) - 1)))
        v_max = std::max(v_max, v.p.dot(loop.lyapunov.X * v.p));
    return v_max;
}

// How far the oracle's verdicts give way to rounding, from the case's own scale. Every trajectory the oracle runs stays
// within rho = sqrt(v_max / lambda_min(X)), so every output inner product along it is at most sigma = c_max rho + h_max.
//  * The construction admits one row on one vertex within an (N + 1)-term inner-product bound on sigma (Higham 2002,
//    Sec. 3.1), and the oracle's evaluation of the same row carries that bound again, so it is counted twice.
//  * The construction's row chain C A_K^k and the oracle's trajectory A_K^k x each add, per product, an N-term bound on
//    ||A_K||_F rho, carried forward by ||A_K^j||_2 <= sqrt(kappa(X)), and in the row chain also through
//    ||A_K^j||_F <= sqrt(N kappa(X)) (Higham 2002, Sec. 3.5). Over the certificate horizon K, with kappa(X) >= 1, both
//    together are within K c_max (1 + sqrt(N)) kappa(X) bound(N, ||A_K||_F rho).
// The level and the certificate are relative to h within a few eps, which the first term covers.
template <std::size_t N>
double depth_of(const closed_loop<N>& loop, double v_max)
{
    auto const& L = loop.lyapunov;
    double const rho = std::sqrt(v_max / L.lambda_min);
    double const c_max = loop.outputs.H.rowwise().norm().maxCoeff();
    double const sigma = c_max * rho + loop.outputs.h.cwiseAbs().maxCoeff();
    double const kappa = L.lambda_max / L.lambda_min;
    auto const horizon = static_cast<double>(std::max<std::int64_t>(certificate_steps(level<N>(loop, 0.0), v_max, L.decay), static_cast<std::int64_t>(N)));
    double const chain = horizon * c_max * (1.0 + std::sqrt(static_cast<double>(N))) * kappa * bound(static_cast<std::int32_t>(N), loop.A_K.norm() * rho);
    return 2.0 * bound(static_cast<std::int32_t>(N + 1), sigma) + chain;
}

// The construction's and the oracle's forward error of one vertex, each carried along the trajectory by sqrt(kappa(X)).
template <std::size_t N>
double vertex_depth(const closed_loop<N>& loop, const vertex<N>& v)
{
    double const c_max = loop.outputs.H.rowwise().norm().maxCoeff();
    return 2.0 * c_max * std::sqrt(loop.lyapunov.lambda_max / loop.lyapunov.lambda_min) * v.error.norm();
}

struct verdict
{
    bool unbounded;
    std::int32_t vertex_outside;
    std::int32_t admissible_outside;
    std::size_t vertices;
    double depth;
};

template <std::size_t N>
bool outside(const rows<N>& set, const point<N>& x)
{
    for(Eigen::Index i = 0; i < set.H.rows(); ++i)
        if(set.H.row(i).dot(x) - set.h(i) > bound(static_cast<std::int32_t>(N + 1), set.H.row(i).cwiseAbs().dot(x.cwiseAbs()) + std::abs(set.h(i))))
            return true;
    return false;
}

// Violations of the returned set S: unbounded, a vertex that leaves Y deepened by the depth, or a point that stays in
// Y made shallower by the depth for all time but lies outside S.
template <std::size_t N>
verdict check(const rows<N>& S, const closed_loop<N>& loop, std::mt19937& g, std::int32_t samples)
{
    verdict result{.unbounded = !positively_spans<N>(S.H), .vertex_outside = 0, .admissible_outside = 0, .vertices = 0, .depth = 0.0};
    if(result.unbounded)
        return result;
    auto const vertices = vertices_of<N>(S);
    auto const box = doubled_box<N>(vertices);
    result.vertices = vertices.size();
    result.depth = depth_of<N>(loop, largest_level<N>(loop, box));
    for(auto const& v : vertices)
        result.vertex_outside += admissible_forever<N>(v.p, loop, result.depth + vertex_depth<N>(loop, v)) ? 0 : 1;
    for(std::int32_t s = 0; s < samples; ++s)
    {
        point<N> x;
        for(Eigen::Index j = 0; j < static_cast<Eigen::Index>(N); ++j)
            x(j) = box.lower(j) + unit_draw(g) * (box.upper(j) - box.lower(j));
        result.admissible_outside += admissible_forever<N>(x, loop, -result.depth) && outside<N>(S, x) ? 1 : 0;
    }
    return result;
}

}

#endif
