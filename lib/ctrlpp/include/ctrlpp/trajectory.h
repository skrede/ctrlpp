#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_H

/// @brief Trajectory generation subsystem -- paths, trajectory segments, and composition.
///
/// Provides elementary polynomial (cubic, quintic, septic) and trigonometric
/// (harmonic, cycloidal) path profiles with variadic piecewise composition.
/// All profiles produce full derivative information (position, velocity, acceleration).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/path_segment.h"
#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/cubic_path.h"
#include "ctrlpp/trajectory/quintic_path.h"
#include "ctrlpp/trajectory/septic_path.h"
#include "ctrlpp/trajectory/harmonic_path.h"
#include "ctrlpp/trajectory/cycloidal_path.h"
#include "ctrlpp/trajectory/trajectory.h"
#include "ctrlpp/trajectory/cubic_trajectory.h"
#include "ctrlpp/trajectory/quintic_trajectory.h"
#include "ctrlpp/trajectory/septic_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/modified_trap_trajectory.h"
#include "ctrlpp/trajectory/modified_sin_trajectory.h"
#include "ctrlpp/trajectory/time_scaling.h"
#include "ctrlpp/trajectory/piecewise_path.h"
#include "ctrlpp/trajectory/piecewise_trajectory.h"
#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/smoothing_spline.h"
#include "ctrlpp/trajectory/bspline_trajectory.h"
#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/synchronize.h"

#endif
