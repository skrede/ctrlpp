#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TINYMPC_SESSION_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_TINYMPC_SESSION_H

#include "json_writer.h"

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <optional>

namespace ctrlpp::bench::qp_solvers
{

/// TinyMPC compiles against the Eigen it bundles, which differs from the one every other unit uses, so
/// nothing of Eigen crosses this boundary: matrices travel as column-major arrays.
struct tinympc_problem
{
    std::int32_t nx;
    std::int32_t nu;
    std::int32_t horizon;
    std::vector<double> transition;
    std::vector<double> actuation;
    std::vector<double> state_weight;
    std::vector<double> input_weight;
    std::vector<double> target;
    double input_limit;
    std::optional<double> state_limit;
};

struct tinympc_options
{
    double tolerance;
    std::int32_t max_iter;
    double rho;
    bool adaptive_rho;
};

/// One TinyMPC solver. TinyMPC's knot count is the horizon plus one: x has N + 1 columns and u has N.
class tinympc_session
{
public:
    tinympc_session(const tinympc_problem& problem, const tinympc_options& options);

    ~tinympc_session();

    tinympc_session(const tinympc_session&) = delete;
    tinympc_session& operator=(const tinympc_session&) = delete;

    void set_initial_state(const double* x0);

    auto solve() -> std::int32_t;

    auto iterations() const -> std::int32_t;

    void read(double* states, double* inputs) const;

private:
    struct solver_state;
    std::unique_ptr<solver_state> state_;
};

auto tinympc_eigen_version() -> std::string;

void describe_tinympc(json_writer& out, const tinympc_options& options);

}

#endif
