#include "PlanarAcousticSolver.h"

::tarch::logging::Log exahype2::training::acoustic::PlanarAcousticSolver::_log(
  "exahype2::training::acoustic::PlanarAcousticSolver"
);

/*
 * Enables the usage of shortcuts to access variables, e.g., use Q[s.p] instead
 * of Q[0]
 */
::exahype2::training::acoustic::VariableShortcuts s;

double ::exahype2::training::acoustic::PlanarAcousticSolver::maxEigenvalue(
  const double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  double                                         dt,
  int                                            normal
) {
  constexpr double K0  = 4.0;
  constexpr double rho = 1.0;
  constexpr double c   = std::sqrt(4.0 / rho);

  return c;
}

void ::exahype2::training::acoustic::PlanarAcousticSolver::flux(
  const double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  double                                         dt,
  int                                            normal,
  double* __restrict__ F
) {
  constexpr double K0  = 4.0;
  constexpr double rho = 1.0;

  switch (normal) {
  case 0: // Flux in x-direction
    F[s.p]     = K0 * Q[s.v + normal];
    F[s.v + 0] = Q[s.p] / rho;
    F[s.v + 1] = 0.0;
    break;
  case 1: // Flux in y-direction
    F[s.p]     = K0 * Q[s.v + normal];
    F[s.v + 0] = 0.0;
    F[s.v + 1] = Q[s.p] / rho;
  }
}

void exahype2::training::acoustic::PlanarAcousticSolver::initialCondition(
  double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  const ::tarch::la::Vector<Dimensions, int>&    index,
  bool                                           gridIsConstructed
) {
  // Simple translation in positive diagonal direction
  const double val = cos(-std::numbers::pi * (x[0] + x[1]));

  Q[s.v + 0] = val;
  Q[s.v + 1] = val;

  constexpr double K0  = 4.0;
  constexpr double rho = 1.0;

  // These are defined by the eigenvector of the plane wave operator
  constexpr double kr = std::sqrt(2 * K0 * rho);
  Q[s.p]              = kr * val;
}

void ::exahype2::training::acoustic::PlanarAcousticSolver::boundaryConditions(
  const double* __restrict__ Qinside,
  double* __restrict__ Qoutside,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  int                                            normal
) {
  // Since we are using periodic boundary conditions, this need never be called.
  assert(false);
}
