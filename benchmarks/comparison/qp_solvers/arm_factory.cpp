#include "osqp_arm.h"
#include "arm_factory.h"
#include "argmin_arms.h"
#include "tinympc_arm.h"
#include "competitor_arms.h"

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto is_argmin(arm_kind kind) -> bool
{
    return kind == arm_kind::argmin_dense_admm || kind == arm_kind::argmin_dense_active_set || kind == arm_kind::argmin_sparse_admm;
}

auto make_argmin_arm(const arm_spec& arm, tier level) -> std::unique_ptr<qp_arm>
{
    if(arm.kind == arm_kind::argmin_dense_active_set)
        return make_argmin_active_set_arm(level);
    if(arm.kind == arm_kind::argmin_dense_admm)
        return make_argmin_dense_admm_arm(level, arm.forced_polish);
    return make_argmin_sparse_admm_arm(level, arm.forced_polish);
}

void describe_argmin_arm(json_writer& out, const arm_spec& arm, tier level)
{
    if(arm.kind == arm_kind::argmin_dense_active_set)
        return describe_argmin_active_set(out, level);
    if(arm.kind == arm_kind::argmin_dense_admm)
        return describe_argmin_dense_admm(out, level, arm.forced_polish);
    describe_argmin_sparse_admm(out, level, arm.forced_polish);
}

}

auto make_arm(const arm_spec& arm, tier level, hpipm_dual_signs hpipm_signs) -> std::unique_ptr<qp_arm>
{
    if(is_argmin(arm.kind))
        return make_argmin_arm(arm, level);
    switch(arm.kind)
    {
        case arm_kind::osqp:
            return make_osqp_arm(level, arm.forced_polish);
        case arm_kind::qpoases:
            return make_qpoases_arm(level);
        case arm_kind::daqp:
            return make_daqp_arm(level);
        case arm_kind::proxqp_dense:
            return make_proxqp_dense_arm(level);
        case arm_kind::proxqp_sparse:
            return make_proxqp_sparse_arm(level);
        case arm_kind::hpipm:
            return make_hpipm_arm(level, hpipm_signs);
        default:
            return make_tinympc_arm(level);
    }
}

void describe_arm(json_writer& out, const arm_spec& arm, tier level)
{
    if(is_argmin(arm.kind))
        return describe_argmin_arm(out, arm, level);
    switch(arm.kind)
    {
        case arm_kind::osqp:
            return describe_osqp(out, level, arm.forced_polish);
        case arm_kind::qpoases:
            return describe_qpoases(out, level);
        case arm_kind::daqp:
            return describe_daqp(out, level);
        case arm_kind::proxqp_dense:
        case arm_kind::proxqp_sparse:
            return describe_proxqp(out, level, arm.kind == arm_kind::proxqp_sparse);
        case arm_kind::hpipm:
            return describe_hpipm(out, level);
        default:
            return describe_tinympc_arm(out, level);
    }
}

}
