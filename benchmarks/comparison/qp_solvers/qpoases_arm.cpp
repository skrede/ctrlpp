#include "competitor_arms.h"

#include <qpOASES.hpp>

#include <optional>
#include <string>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

using row_major = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

auto qpoases_options(tier level) -> qpOASES::Options
{
    qpOASES::Options options;
    if(level == tier::mpc)
        options.setToMPC();
    else
        options.setToReliable();
    options.printLevel = qpOASES::PL_NONE;
    return options;
}

auto status_token(qpOASES::returnValue status) -> std::string
{
    switch(status)
    {
        case qpOASES::SUCCESSFUL_RETURN:
            return "SUCCESSFUL_RETURN";
        case qpOASES::RET_MAX_NWSR_REACHED:
            return "RET_MAX_NWSR_REACHED";
        case qpOASES::RET_INIT_FAILED_INFEASIBILITY:
            return "RET_INIT_FAILED_INFEASIBILITY";
        case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY:
            return "RET_HOTSTART_STOPPED_INFEASIBILITY";
        default:
            return "RET_" + std::to_string(static_cast<int>(status));
    }
}

auto claim_of(qpOASES::returnValue status) -> solver_claim
{
    switch(qpOASES::getSimpleStatus(status))
    {
        case 0:
            return solver_claim::solved;
        case 1:
            return solver_claim::iteration_cap;
        case -2:
        case -3:
            return solver_claim::infeasible;
        default:
            return solver_claim::error;
    }
}

/// qpOASES reports a multiplier as positive on an active lower bound, the opposite of OSQP's sign;
/// its bound multipliers come first, then its constraint multipliers. It keeps pointers to the
/// matrices it was initialized with rather than copies, so the arm owns them.
class qpoases_arm final : public qp_arm
{
public:
    explicit qpoases_arm(const qpOASES::Options& options)
        : H_{}
        , A_{}
        , g_{}
        , lb_{}
        , ub_{}
        , boxes_(0)
        , lbA_{}
        , ubA_{}
        , dual_{}
        , options_(options)
        , work_sets_(0)
        , status_(qpOASES::SUCCESSFUL_RETURN)
        , solver_{}
    {
    }

    void setup(const qp_instance& instance, const qp_tick& tick) override
    {
        const auto n = static_cast<qpOASES::int_t>(instance.P.rows());
        boxes_ = instance.box_rows;
        const auto constraints = static_cast<qpOASES::int_t>(instance.A.rows() - boxes_);
        H_ = instance.P;
        A_ = instance.A.bottomRows(constraints);
        solver_.emplace(n, constraints);
        solver_->setOptions(options_);
        lb_ = Eigen::VectorXd::Constant(n, -qpOASES::INFTY);
        ub_ = Eigen::VectorXd::Constant(n, qpOASES::INFTY);
        update(tick);
        work_sets_ = iteration_cap;
        status_ = solver_->init(H_.data(), g_.data(), A_.data(), bounds(lb_), bounds(ub_), lbA_.data(), ubA_.data(), work_sets_);
    }

    void update(const qp_tick& tick) override
    {
        g_ = tick.q;
        lb_.head(boxes_) = tick.l.head(boxes_);
        ub_.head(boxes_) = tick.u.head(boxes_);
        lbA_ = tick.l.tail(tick.l.size() - boxes_);
        ubA_ = tick.u.tail(tick.u.size() - boxes_);
    }

    void solve() override
    {
        work_sets_ = iteration_cap;
        status_ = solver_->hotstart(g_.data(), bounds(lb_), bounds(ub_), lbA_.data(), ubA_.data(), work_sets_);
    }

    void answer(arm_answer& out) override;

private:
    row_major H_;
    row_major A_;
    Eigen::VectorXd g_;
    Eigen::VectorXd lb_;
    Eigen::VectorXd ub_;
    std::int32_t boxes_;
    Eigen::VectorXd lbA_;
    Eigen::VectorXd ubA_;
    Eigen::VectorXd dual_;
    qpOASES::Options options_;
    qpOASES::int_t work_sets_;
    qpOASES::returnValue status_;
    std::optional<qpOASES::QProblem> solver_;

    /// A problem without box rows is posed with no variable bounds at all.
    auto bounds(const Eigen::VectorXd& values) const -> const double*
    {
        return boxes_ > 0 ? values.data() : nullptr;
    }
};

void qpoases_arm::answer(arm_answer& out)
{
    const Eigen::Index n = g_.size();
    const Eigen::Index constraints = lbA_.size();
    out.x.resize(n);
    dual_.resize(n + constraints);
    solver_->getPrimalSolution(out.x.data());
    solver_->getDualSolution(dual_.data());
    out.y.resize(boxes_ + constraints);
    out.y.head(boxes_) = -dual_.head(boxes_);
    out.y.tail(constraints) = -dual_.tail(constraints);
    out.iterations = static_cast<std::int32_t>(work_sets_);
    out.status = status_token(status_);
    out.claim = claim_of(status_);
    out.fallback = false;
}

}

auto make_qpoases_arm(tier level) -> std::unique_ptr<qp_arm>
{
    return std::make_unique<qpoases_arm>(qpoases_options(level));
}

void describe_qpoases(json_writer& out, tier level)
{
    const qpOASES::Options o = qpoases_options(level);
    out.text("preset", level == tier::mpc ? "MPC" : "RELIABLE");
    out.text("print_level", o.printLevel == qpOASES::PL_NONE ? "NONE" : "OTHER");
    out.integer("nwsr_cap", iteration_cap);
    out.flag("enableRamping", o.enableRamping == qpOASES::BT_TRUE);
    out.flag("enableFarBounds", o.enableFarBounds == qpOASES::BT_TRUE);
    out.flag("enableFlippingBounds", o.enableFlippingBounds == qpOASES::BT_TRUE);
    out.flag("enableRegularisation", o.enableRegularisation == qpOASES::BT_TRUE);
    out.flag("enableFullLITests", o.enableFullLITests == qpOASES::BT_TRUE);
    out.flag("enableNZCTests", o.enableNZCTests == qpOASES::BT_TRUE);
    out.integer("enableDriftCorrection", o.enableDriftCorrection);
    out.integer("enableCholeskyRefactorisation", o.enableCholeskyRefactorisation);
    out.flag("enableEqualities", o.enableEqualities == qpOASES::BT_TRUE);
    out.number("terminationTolerance", o.terminationTolerance);
    out.number("boundTolerance", o.boundTolerance);
    out.number("boundRelaxation", o.boundRelaxation);
    out.number("epsNum", o.epsNum);
    out.number("epsDen", o.epsDen);
    out.number("epsRegularisation", o.epsRegularisation);
    out.integer("numRefinementSteps", o.numRefinementSteps);
}

}
