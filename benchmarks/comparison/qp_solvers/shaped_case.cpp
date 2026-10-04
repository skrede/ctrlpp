#include "cases.h"

#include "qp/qp_accuracy.h"

namespace ctrlpp::bench::qp_solvers
{

namespace
{

class shaped_case final : public problem_case
{
public:
    shaped_case()
        : program_(problems::qp::make_dense_mpc_program())
        , instance_{program_.P, program_.A, program_.P.sparseView(), program_.A.sparseView(), 0, std::nullopt}
    {
    }

    auto instance() const -> const qp_instance& override
    {
        return instance_;
    }

    void fill(const Eigen::VectorXd&, qp_tick& tick) const override
    {
        tick.q = program_.q;
        tick.l = program_.l;
        tick.u = program_.u;
        tick.x0.resize(0);
    }

    auto rows_of(const qp_tick& tick, const Eigen::VectorXd& x) const -> inequality_rows override
    {
        return {program_.A * x, tick.l, tick.u};
    }

    auto state_dimension() const -> std::int32_t override
    {
        return 0;
    }

    auto amplitude() const -> double override
    {
        return 0.0;
    }

private:
    problems::qp::dense_program program_;
    qp_instance instance_;
};

}

auto make_shaped_case() -> std::unique_ptr<problem_case>
{
    return std::make_unique<shaped_case>();
}

}
