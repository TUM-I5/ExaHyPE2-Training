#include "EulerSolver.h"

::tarch::logging::Log exahype2::training::euler::EulerSolver::_log(
  "exahype2::training::euler::EulerSolver"
);

::exahype2::training::euler::VariableShortcuts s;

// c chord length, t maximum relative thickness.
double airfoilSymmetric(double x, double c = 100, double t = 0.12) {
  return 5.0 * t * c
         * (0.2969 * sqrt(x / c) + ((((-0.1015) * (x / c) + 0.2843) * (x / c) - 0.3516) * (x / c) - 0.1260) * (x / c));
}

double airfoilCambered(
  double x,
  double c = 100,
  double p = 0.6,
  double m = 0.04,
  double t = 0.12
) {
  double yc; // y position along camber line
  double divisor;
  if (0.0 <= x / c && x / c <= p) {
    yc      = m * x * (2.0 * p - x / c) / p / p;
    divisor = p * p;
  } else if (p <= x / c && x / c <= 1.0) {
    yc      = m * (x - c) * (2.0 * p - x / c - 1.0) / (1.0 - p) / (1.0 - p);
    divisor = pow(1.0 - p, 2);
  } else {
    yc      = 0.0;
    divisor = 1.0;
  }

  const double dycdx = 2.0 * m * (p - x / c) / divisor;
  const double theta = atan(dycdx);

  const double yt
    = 5.0 * t * c
      * (0.2969 * sqrt(x / c) + ((((-0.1015) * (x / c) + 0.2843) * (x / c) - 0.3516) * (x / c) - 0.1260) * (x / c));

  double xu = x - yt * sin(theta);
  double yu = yc + yt * cos(theta);
  double xl = x + yt * sin(theta);
  double yl = yc - yt * cos(theta);

  return 1.0;
}

void exahype2::training::euler::EulerSolver::initialCondition(
  double* __restrict__ Q,
  const tarch::la::Vector<Dimensions, double>& x,
  const tarch::la::Vector<Dimensions, double>& h,
  bool                                         gridIsConstructed
) {
  const double y     = airfoilSymmetric(x[0]);
  const double alpha = (y > x[1] && -y < x[1]) ? 0.0 : 1.0;

  Q[0] = (y > x[1] && -y < x[1]) ? 1000.0 : 1.0;
  Q[1] = alpha * 1.0;
  Q[2] = 0.0;
  Q[3] = alpha * 2.5;
}

void exahype2::training::euler::EulerSolver::boundaryConditions(
  const double* __restrict__ Qinside,
  double* __restrict__ Qoutside,
  const tarch::la::Vector<Dimensions, double>& x,
  const tarch::la::Vector<Dimensions, double>& h,
  double                                       t,
  int                                          normal
) {
  if (normal == 0 && x[0] == DomainOffset[0]) {
    Qoutside[0] = 1.0;
    Qoutside[1] = 1.0;
    Qoutside[2] = 0.0;
    Qoutside[3] = 2.5;
  } else {
    Qoutside[0] = Qinside[0];
    Qoutside[1] = Qinside[1];
    Qoutside[2] = Qinside[2];
    Qoutside[3] = Qinside[3];
  }
}

double ::exahype2::training::euler::EulerSolver::maxEigenvalue(
  const double* __restrict__ Q,
  const tarch::la::Vector<Dimensions, double>& x,
  const tarch::la::Vector<Dimensions, double>& h,
  double                                       t,
  double                                       dt,
  int                                          normal
) {
  const double     irho  = 1.0 / Q[0];
  constexpr double gamma = 1.4;
  const double     p     = (gamma - 1)
                   * (Q[3] - 0.5 * irho * (Q[1] * Q[1] + Q[2] * Q[2]));

  const double c = std::sqrt(gamma * p * irho);
  const double u = Q[normal + 1] * irho;

  return std::max(std::abs(u - c), std::abs(u + c));
}

void ::exahype2::training::euler::EulerSolver::flux(
  const double* __restrict__ Q,
  const tarch::la::Vector<Dimensions, double>& x,
  const tarch::la::Vector<Dimensions, double>& h,
  double                                       t,
  double                                       dt,
  int                                          normal,
  double* __restrict__ F
) {
  const double     irho  = 1.0 / Q[0];
  constexpr double gamma = 1.4;
  const double     p     = (gamma - 1)
                   * (Q[3] - 0.5 * irho * (Q[1] * Q[1] + Q[2] * Q[2]));

  F[0] = Q[normal + 1];
  F[1] = Q[normal + 1] * Q[1] * irho;
  F[2] = Q[normal + 1] * Q[2] * irho;
  F[3] = Q[normal + 1] * irho * (Q[3] + p);

  F[normal + 1] += p;
}
