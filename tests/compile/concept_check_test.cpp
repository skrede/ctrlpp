#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/particle_filter.h"
#include "ctrlpp/estimation/observer_policy.h"
#include "ctrlpp/model/differentiable_dynamics.h"
#include "ctrlpp/model/differentiable_measurement.h"
#include "ctrlpp/model/measurement_model.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/estimation/sigma_points/sigma_point_strategy.h"
#include "ctrlpp/estimation/sigma_points/merwe_sigma_points.h"
#include "ctrlpp/estimation/sigma_points/julier_sigma_points.h"
#include "ctrlpp/estimation/resampling/resampling_strategy.h"
#include "ctrlpp/estimation/resampling/systematic_resampling.h"
#include "ctrlpp/estimation/resampling/multinomial_resampling.h"

#include <random>

using namespace ctrlpp;

// -- measurement_model --

// A lambda-like callable satisfying measurement_model<M, double, 2, 1>
struct LinearMeasurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1> {
        return Vector<double, 1>{x[0]};
    }
};

static_assert(measurement_model<LinearMeasurement, double, 2, 1>);

// Linear C matrix wrapper trivially satisfies measurement_model
struct CMatrixWrapper {
    Matrix<double, 1, 2> C{Matrix<double, 1, 2>::Identity()};

    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1> {
        return C * x;
    }
};

static_assert(measurement_model<CMatrixWrapper, double, 2, 1>);

// -- dynamics_model --

struct SimpleDynamics {
    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2> {
        Vector<double, 2> result;
        result[0] = x[0] + u[0];
        result[1] = x[1];
        return result;
    }
};

static_assert(dynamics_model<SimpleDynamics, double, 2, 1>);

// -- differentiable_dynamics --

struct DiffDynamics {
    auto operator()(const Vector<double, 2>& x,
                    const Vector<double, 1>& u) const -> Vector<double, 2> {
        Vector<double, 2> result;
        result[0] = x[0] + u[0];
        result[1] = x[1];
        return result;
    }

    auto jacobian_x(const Vector<double, 2>&,
                    const Vector<double, 1>&) const -> Matrix<double, 2, 2> {
        return Matrix<double, 2, 2>::Identity();
    }

    auto jacobian_u(const Vector<double, 2>&,
                    const Vector<double, 1>&) const -> Matrix<double, 2, 1> {
        Matrix<double, 2, 1> B;
        B << 1.0, 0.0;
        return B;
    }
};

static_assert(differentiable_dynamics<DiffDynamics, double, 2, 1>);

// SimpleDynamics satisfies dynamics_model but NOT differentiable_dynamics
static_assert(dynamics_model<SimpleDynamics, double, 2, 1>);
static_assert(!differentiable_dynamics<SimpleDynamics, double, 2, 1>);

// -- differentiable_measurement --

struct DiffMeasurement {
    auto operator()(const Vector<double, 2>& x) const -> Vector<double, 1> {
        return Vector<double, 1>{x[0]};
    }

    auto jacobian(const Vector<double, 2>&) const -> Matrix<double, 1, 2> {
        Matrix<double, 1, 2> H;
        H << 1.0, 0.0;
        return H;
    }
};

static_assert(differentiable_measurement<DiffMeasurement, double, 2, 1>);

// LinearMeasurement satisfies measurement_model but NOT differentiable_measurement
static_assert(measurement_model<LinearMeasurement, double, 2, 1>);
static_assert(!differentiable_measurement<LinearMeasurement, double, 2, 1>);

// -- UKF observer concept checks --

using ukf_test_t = ukf<double, 2, 1, 1, SimpleDynamics, LinearMeasurement>;
static_assert(ObserverPolicy<ukf_test_t>);
static_assert(CovarianceObserver<ukf_test_t>);

// -- Particle filter observer concept checks --

using pf_test_t = particle_filter<double, 2, 1, 1, 10, SimpleDynamics, LinearMeasurement>;
static_assert(ObserverPolicy<pf_test_t>);
static_assert(!CovarianceObserver<pf_test_t>);

// -- Sigma point strategy concept checks --

static_assert(sigma_point_strategy<merwe_sigma_points<double, 2>, double, 2>);
static_assert(sigma_point_strategy<julier_sigma_points<double, 2>, double, 2>);

// -- Resampling strategy concept checks --

static_assert(resampling_strategy<systematic_resampling, std::mt19937_64, 10>);
static_assert(resampling_strategy<multinomial_resampling, std::mt19937_64, 10>);

int main() {
    return 0;
}
