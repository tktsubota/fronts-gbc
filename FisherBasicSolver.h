#include <armadillo>
#include "InitialCondition.h"
#include "Coefficient.h"
#include <iostream>
#include <H5Cpp.h>
#include "H5DatasetBatchLoader.h"

#ifndef FISHERBASICSOLVER_H
#define FISHERBASICSOLVER_H

/* u_t = a(t)*u(1-u) + d(t)u_xx + c(t)u_x = a(t)u + d(t)u_xx + c(t)u_x - a(t)u^2 */
class FisherBasicSolver {
    protected:
        Coefficient *a;
        Coefficient *d;
        Coefficient *c;
        double t0, t1, dt;
        double x0, x1, dx;
        InitialCondition *ic;
        size_t tout;
        double buf; // boundary condition buffer region size (naive bc if zero)
        double gamma; // frame speed multiplier (c = gamma * d)
        std::string impulse_file; // file with impulse
        std::string ic_integral_file; // file with the IC integral

        size_t Nt;
        size_t Nx;

        arma::vec t_full; // all evaluated t coords
        arma::vec t; // output t coords (considering tout)
        arma::vec x; // output x coords

        bool use_td_bc;
        bool use_dirichlet; // False uses naive zero Neumann regardless of use_td_bc
        bool use_inline_greens;     // compute I0/I1/I2 on-the-fly (no file)
        bool use_inline_ic_integral; // compute IC integral on-the-fly (no file)

        std::unique_ptr<H5DatasetBatchLoader> I0_loader;
        std::unique_ptr<H5DatasetBatchLoader> I1_loader;
        std::unique_ptr<H5DatasetBatchLoader> I2_loader;
        arma::vec a_all;
        arma::vec d_all;
        arma::vec c_all;

        arma::vec IC_integral;

        // Construct matrix for implicit numerical scheme
        arma::sp_mat M;
        void init_matrix();
        void update_matrix(double t_curr);

        // Define right hand side solution vector
        arma::vec b;
        virtual void update_rhs(arma::vec u, double t_curr);

        // Compute boundary value
        double compute_boundary_value(double t, size_t k, arma::vec g, arma::vec gp);

        // Precompute I0/I1/I2 matrices to file, then open via batch loaders
        void precompute_greens_to_file();
        // Precompute IC integral to file, then load into IC_integral
        void precompute_ic_integral_to_file();

    public:
        struct Result {
            arma::vec t;
            arma::vec x;
            arma::mat u;
        };

        FisherBasicSolver(Coefficient *a, Coefficient *d, Coefficient *c,
                          double t0, double t1, double dt,
                          double x0, double x1, double dx,
                          InitialCondition *ic, size_t tout,
                          double buf, double gamma,
                          std::string impulse_file, std::string ic_integral_file);

        Result run(std::ostream &logger = std::cout, bool verbose = true);

        virtual ~FisherBasicSolver() = default;

    };

#endif // FISHERBASICSOLVER_H
