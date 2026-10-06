#include "Fisher23Solver.h"

void Fisher23Solver::update_rhs(arma::vec u, double t_curr) {
    // Neumann boundary condition on left
    b(0) = 0;

    for (size_t i=1; i < Nx-1; i += 1){
        b(i) = u(i) + ((u(i) * u(i)) - (1 + a->func(t_curr))*(u(i) * u(i) * u(i)))*dt;
    }
}
