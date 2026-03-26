#ifndef HPP_GUARD_CTRLPP_TRAJ_H
#define HPP_GUARD_CTRLPP_TRAJ_H

/// @brief Trajectory generation subsystem -- paths, trajectory segments, and composition.
///
/// Provides elementary polynomial (cubic, quintic, septic) and trigonometric
/// (harmonic, cycloidal) path profiles with variadic piecewise composition.
/// All profiles produce full derivative information (position, velocity, acceleration).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009

#include "ctrlpp/traj/trajectory_types.h"
#include "ctrlpp/traj/path_segment.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/cubic_path.h"
#include "ctrlpp/traj/quintic_path.h"
#include "ctrlpp/traj/septic_path.h"
#include "ctrlpp/traj/harmonic_path.h"
#include "ctrlpp/traj/cycloidal_path.h"
#include "ctrlpp/traj/trajectory.h"
#include "ctrlpp/traj/cubic_trajectory.h"
#include "ctrlpp/traj/quintic_trajectory.h"
#include "ctrlpp/traj/septic_trajectory.h"
#include "ctrlpp/traj/modified_trap_trajectory.h"
#include "ctrlpp/traj/piecewise_path.h"
#include "ctrlpp/traj/piecewise_trajectory.h"

#endif
