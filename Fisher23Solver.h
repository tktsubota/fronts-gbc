#pragma once

#include "FisherBasicSolver.h"

#ifndef FISHER23SOLVER_H
#define FISHER23SOLVER_H

/* u_t = a(t)*u + d(t)u_xx + c(t)u_x + u^2 - (1 + a(t))u^3 */
class Fisher23Solver: public FisherBasicSolver {
public:
    using FisherBasicSolver::FisherBasicSolver; // inherit constructors

protected:
    void update_rhs(arma::vec u, double t_curr) override;
};

#endif // FISHER23SOLVER_H
