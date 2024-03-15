#include "SWESolver.h"

::tarch::logging::Log exahype2::training::swe::SWESolver::_log(
  "exahype2::aderdg::swe::SWESolver"
);

/*
 * Enables the usage of shortcuts to access variables, e.g., use Q[s.h] instead
 * of Q[0]
 */
::exahype2::training::swe::VariableShortcuts s;

void ::exahype2::training::swe::SWESolver::initialCondition(
  double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  const ::tarch::la::Vector<Dimensions, int>&    index,
  bool                                           gridIsConstructed
) {
  Q[s.hu + 0] = 0.0; // v_x
  Q[s.hu + 1] = 0.0; // v_y

  // Part 1
  /*
  Q[s.h] = x[0] < 0.0 ? 1.0 : 2.0; // h
  Q[s.b] = 0.0;
  */

  // Part 2
  Q[s.h] = 1.0;
  Q[s.b] = (tarch::la::norm2(x) < 0.5 ? 0.2 : 0.0);
}

void ::exahype2::training::swe::SWESolver::boundaryConditions(
  const double* __restrict__ Qinside,
  double* __restrict__ Qoutside,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  int                                            normal
) {
  Qoutside[s.h]      = 1.0; // h
  Qoutside[s.hu + 0] = 0.0; // v_x
  Qoutside[s.hu + 1] = 0.0; // v_y
  Qoutside[s.b]      = 0.0; // b
}

double ::exahype2::training::swe::SWESolver::maxEigenvalue(
  const double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  double                                         dt,
  int                                            normal
) {
  constexpr double grav = 9.81;

  const double u = Q[s.hu + normal] / Q[s.h];
  const double c = std::sqrt(grav * Q[s.h]);

  return std::max(std::abs(u + c), std::abs(u - c));
}

void ::exahype2::training::swe::SWESolver::flux(
  const double* __restrict__ Q,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  double                                         dt,
  int                                            normal,
  double* __restrict__ F
) {
  double ih = 1.0 / Q[s.h];

  // Part 1
  /*
  constexpr double grav = 9.81;
  F[s.h]      = Q[s.hu + normal];
  F[s.hu + 0] = Q[s.hu + normal] * Q[s.hu + 0] * ih;
  F[s.hu + 1] = Q[s.hu + normal] * Q[s.hu + 1] * ih;
  F[s.b]      = 0.0;
  F[s.hu + normal] += 0.5 * grav * Q[s.h] * Q[s.h];
  */

  // Part 2
  F[s.h]      = Q[s.hu + normal];
  F[s.hu + 0] = Q[s.hu + normal] * Q[s.hu + 0] * ih;
  F[s.hu + 1] = Q[s.hu + normal] * Q[s.hu + 1] * ih;
  F[s.b]      = 0.0;
}

void ::exahype2::training::swe::SWESolver::nonconservativeProduct(
  const double* __restrict__ Q,
  const double* __restrict__ deltaQ,
  const ::tarch::la::Vector<Dimensions, double>& x,
  const ::tarch::la::Vector<Dimensions, double>& h,
  double                                         t,
  double                                         dt,
  int                                            normal,
  double* __restrict__ BTimesDeltaQ
) {
  // Part 1
  /*
  BTimesDeltaQ[s.h]      = 0.0;
  BTimesDeltaQ[s.hu + 0] = 0.0;
  BTimesDeltaQ[s.hu + 1] = 0.0;
  BTimesDeltaQ[s.b]      = 0.0;
  */

  // Part 2
  constexpr double grav = 9.81;
  BTimesDeltaQ[s.h] = 0.0;
  switch (normal) {
  case 0:
    BTimesDeltaQ[s.hu + 0] = grav * Q[s.h] * (deltaQ[s.h] + deltaQ[s.b]);
    BTimesDeltaQ[s.hu + 1] = 0.0;
    break;
  case 1:
    BTimesDeltaQ[s.hu + 0] = 0.0;
    BTimesDeltaQ[s.hu + 1] = grav * Q[s.h] * (deltaQ[s.h] + deltaQ[s.b]);
    break;
  }
  BTimesDeltaQ[s.b] = 0.0;
}

::exahype2::RefinementCommand exahype2::training::swe::SWESolver::
  refinementCriterion(
    const double* __restrict__ Q,
    const ::tarch::la::Vector<Dimensions, double>& x,
    const ::tarch::la::Vector<Dimensions, double>& h,
    double                                         t
  ) {
  // Part 1
  //return ::exahype2::RefinementCommand::Keep;

  // Part 2
  return (std::abs(tarch::la::norm2(x) - 0.5) < 0.2 ?
    ::exahype2::RefinementCommand::Refine :
    ::exahype2::RefinementCommand::Keep);
}
