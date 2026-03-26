#ifndef HPP_GUARD_CTRLPP_TRAJ_H
#define HPP_GUARD_CTRLPP_TRAJ_H

/// @brief Trajectory generation subsystem -- motion laws, segments, and composition.
///
/// Provides elementary polynomial (cubic, quintic, septic) and trigonometric
/// (harmonic, cycloidal) trajectory profiles with variadic piecewise composition.
/// All profiles produce full derivative information (position, velocity, acceleration).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009

#include "ctrlpp/traj/trajectory_types.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/motion_laws.h"
#include "ctrlpp/traj/scaled_segment.h"
#include "ctrlpp/traj/cubic_segment.h"
#include "ctrlpp/traj/quintic_segment.h"
#include "ctrlpp/traj/septic_segment.h"
#include "ctrlpp/traj/piecewise.h"

#endif
