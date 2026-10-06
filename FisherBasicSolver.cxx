#include "FisherBasicSolver.h"
#include "GreensIntegrals.h"
#include "options.h"
#include <fstream>
#include <vector>

static bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

// -----------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------

FisherBasicSolver::FisherBasicSolver(
        Coefficient *a, Coefficient *d, Coefficient *c,
        double t0, double t1, double dt,
        double x0, double x1, double dx,
        InitialCondition *ic, size_t tout,
        double buf, double gamma,
        std::string impulse_file, std::string ic_integral_file)
    : a(a), d(d), c(c),
      t0(t0), t1(t1), dt(dt),
      x0(x0), x1(x1), dx(dx),
      ic(ic), tout(tout),
      buf(buf), gamma(gamma),
      impulse_file(impulse_file),
      ic_integral_file(ic_integral_file),
      use_td_bc(false), use_dirichlet(true),
      use_inline_greens(false), use_inline_ic_integral(false)
{
    Nt = (size_t) round((t1 - t0) / (dt*tout));
    Nx = (size_t) round((x1 - x0) / dx) + 1;

    t_full = arma::linspace(t0, t1-dt, Nt*tout);
    t      = arma::linspace(t0, t1-(dt*tout), Nt);
    x      = arma::linspace(x0, x1, Nx);

    if (buf <= 0) {
        if (buf <= -1) use_dirichlet = false;
        return;
    }

    use_td_bc = true;

    // Compute a_all, d_all, c_all (needed regardless of mode)
    a_all = arma::zeros<arma::vec>(Nt*tout);
    d_all = arma::zeros<arma::vec>(Nt*tout);
    c_all = arma::zeros<arma::vec>(Nt*tout);
    for (size_t i = 0; i < Nt*tout; i++) {
        a_all[i] = a->func(t_full[i]);
        d_all[i] = d->func(t_full[i]);
        c_all[i] = c->func(t_full[i]);
    }

    // ---- I0/I1/I2: inline, load, or precompute ----
    if (impulse_file.empty()) {
        use_inline_greens = true;
        std::cout << "Using inline Green's function computation." << std::endl;
    } else if (file_exists(impulse_file)) {
        H5::H5File file(impulse_file, H5F_ACC_RDONLY);
        size_t batch_size = 100;
        I0_loader = std::make_unique<H5DatasetBatchLoader>(file, "I0", batch_size);
        I1_loader = std::make_unique<H5DatasetBatchLoader>(file, "I1", batch_size);
        I2_loader = std::make_unique<H5DatasetBatchLoader>(file, "I2", batch_size);
        std::cout << "Loaded impulse file " << impulse_file << std::endl;
    } else {
        std::cout << "Precomputing Green's function integrals to " << impulse_file << " ..." << std::endl;
        precompute_greens_to_file();
        H5::H5File file(impulse_file, H5F_ACC_RDONLY);
        size_t batch_size = 100;
        I0_loader = std::make_unique<H5DatasetBatchLoader>(file, "I0", batch_size);
        I1_loader = std::make_unique<H5DatasetBatchLoader>(file, "I1", batch_size);
        I2_loader = std::make_unique<H5DatasetBatchLoader>(file, "I2", batch_size);
        std::cout << "Loaded precomputed impulse file " << impulse_file << std::endl;
    }

    // ---- IC integral: inline, load, or precompute ----
    if (ic_integral_file.empty()) {
        use_inline_ic_integral = true;
        std::cout << "Using inline IC integral computation." << std::endl;
    } else if (file_exists(ic_integral_file)) {
        H5::H5File file_ic(ic_integral_file, H5F_ACC_RDONLY);
        H5::DataSet dataset = file_ic.openDataSet("I_ic");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims);
        std::vector<double> buffer(dims[0]);
        dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
        IC_integral = arma::vec(buffer.data(), dims[0], true);
        std::cout << "Loaded ic file " << ic_integral_file << std::endl;
    } else {
        std::cout << "Precomputing IC integral to " << ic_integral_file << " ..." << std::endl;
        precompute_ic_integral_to_file();
        H5::H5File file_ic(ic_integral_file, H5F_ACC_RDONLY);
        H5::DataSet dataset = file_ic.openDataSet("I_ic");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims);
        std::vector<double> buffer(dims[0]);
        dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
        IC_integral = arma::vec(buffer.data(), dims[0], true);
        std::cout << "Loaded precomputed ic file " << ic_integral_file << std::endl;
    }
}

// -----------------------------------------------------------------------
// Precompute helpers
// -----------------------------------------------------------------------

void FisherBasicSolver::precompute_greens_to_file() {
    size_t N = Nt * tout;
    // Allocate three N×N matrices (row i, col j = I0(buf, t[i], t[j]))
    std::vector<double> I0_data(N * N, 0.0);
    std::vector<double> I1_data(N * N, 0.0);
    std::vector<double> I2_data(N * N, 0.0);

    size_t report_every = N / 20;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            // row-major storage: [i*N + j]
            I0_data[i * N + j] = I0_inline(buf, t_full[i], t_full[j], gamma, *a, *d);
            I1_data[i * N + j] = I1_inline(buf, t_full[i], t_full[j], gamma, *a, *d);
            I2_data[i * N + j] = I2_inline(buf, t_full[i], t_full[j], gamma, *a, *d);
        }
        if (report_every > 0 && (i+1) % report_every == 0) {
            std::cout << "  Precompute progress: " << (i+1)*100/N << "%" << std::endl;
        }
    }

    H5::H5File file(impulse_file, H5F_ACC_TRUNC);

    hsize_t dims2[2] = {(hsize_t)N, (hsize_t)N};
    H5::DataSpace ds2(2, dims2);
    file.createDataSet("I0", H5::PredType::NATIVE_DOUBLE, ds2).write(I0_data.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("I1", H5::PredType::NATIVE_DOUBLE, ds2).write(I1_data.data(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("I2", H5::PredType::NATIVE_DOUBLE, ds2).write(I2_data.data(), H5::PredType::NATIVE_DOUBLE);

    hsize_t dims1[1] = {(hsize_t)N};
    H5::DataSpace ds1(1, dims1);
    file.createDataSet("t", H5::PredType::NATIVE_DOUBLE, ds1).write(t_full.memptr(), H5::PredType::NATIVE_DOUBLE);

    file.createAttribute("buffer_size", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &buf);
    file.createAttribute("dt",          H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &dt);
    file.createAttribute("T",           H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &t1);
    file.close();
}

void FisherBasicSolver::precompute_ic_integral_to_file() {
    size_t N = Nt * tout;
    std::vector<double> vals(N, 0.0);

    ExponentialIC* eic = dynamic_cast<ExponentialIC*>(ic);
    if (eic) {
        double x_boundary = x1 - buf;
        for (size_t k = 0; k < N; k++) {
            vals[k] = ic_int_inline(buf, t_full[k], gamma, *a, *d, x_boundary, eic->lambda);
        }
    }

    H5::H5File file(ic_integral_file, H5F_ACC_TRUNC);

    hsize_t dims1[1] = {(hsize_t)N};
    H5::DataSpace ds1(1, dims1);
    file.createDataSet("t",    H5::PredType::NATIVE_DOUBLE, ds1).write(t_full.memptr(), H5::PredType::NATIVE_DOUBLE);
    file.createDataSet("I_ic", H5::PredType::NATIVE_DOUBLE, ds1).write(vals.data(),     H5::PredType::NATIVE_DOUBLE);

    file.createAttribute("buffer_size", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &buf);
    file.createAttribute("dt",          H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &dt);
    file.createAttribute("T",           H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &t1);
    file.createAttribute("x1",          H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &x1);
    file.close();
}

// -----------------------------------------------------------------------
// Matrix / RHS setup (unchanged)
// -----------------------------------------------------------------------

void FisherBasicSolver::init_matrix() {
    const size_t nnz = 3 * Nx - 3;

    arma::umat locs(2, nnz);
    arma::vec  vals(nnz, arma::fill::zeros);

    size_t k = 0;

    for (size_t i = 1; i < Nx - 1; ++i) {
        locs(0, k) = i; locs(1, k) = i - 1; k++;
        locs(0, k) = i; locs(1, k) = i;     k++;
        locs(0, k) = i; locs(1, k) = i + 1; k++;
    }

    locs(0, k) = 0; locs(1, k) = 0; k++;
    locs(0, k) = 0; locs(1, k) = 1; k++;

    locs(0, k) = Nx - 1; locs(1, k) = Nx - 1; k++;

    M = arma::sp_mat(locs, vals, Nx, Nx);
}

void FisherBasicSolver::update_matrix(double t_curr) {
    double diff = d->func(t_curr) * dt / (dx*dx);
    double adv  = (c->func(t_curr) * dt) / (2.0*dx);
    double diag_val = (1.0 - a->func(t_curr) * dt);

    for (size_t i = 1; i < Nx - 1; ++i) {
        M(i, i-1) = -diff + adv;
        M(i, i)   = diag_val + 2.0*diff;
        M(i, i+1) = -diff - adv;
    }

    M(0,0) = 1.0;
    M(0,1) = -1.0;

    if (use_dirichlet) {
        M(Nx-1, Nx-1) = 1.0;
    } else {
        M(Nx-1, Nx-2) = 1.0;
        M(Nx-1, Nx-1) = -1.0;
    }
}

void FisherBasicSolver::update_rhs(arma::vec u, double t_curr) {
    b(0) = 0;

    for (size_t i=1; i < Nx-1; i += 1){
        b(i) = u(i) - a->func(t_curr)*(u(i) * u(i))*dt;
    }
}

// -----------------------------------------------------------------------
// Boundary value computation
// -----------------------------------------------------------------------

double trapezoidal_integrate(arma::vec x, arma::vec y, double x0, double x1, double dx) {
    double ans = 0;
    size_t idx = 0;
    while (x(idx) < x0) {
        idx += 1;
    }
    ans += y(idx);
    idx += 1;
    for (; x(idx) < x1; idx++) {
        ans += 2*y(idx);
    }
    ans += y(idx);
    ans *= dx/2;
    return ans;
}

double FisherBasicSolver::compute_boundary_value(double t_curr, size_t k, arma::vec g, arma::vec gp) {
    arma::vec h0 = (a_all - 2*d_all) % g - gp;
    arma::vec h1 = -2*c_all % g;
    arma::vec h2 = 4*d_all % g;

    arma::vec col0, col1, col2;

    if (use_inline_greens) {
        col0 = arma::zeros<arma::vec>(t_full.n_elem);
        col1 = arma::zeros<arma::vec>(t_full.n_elem);
        col2 = arma::zeros<arma::vec>(t_full.n_elem);
        for (size_t j = 0; j < k; j++) {
            col0[j] = I0_inline(buf, t_full[k], t_full[j], gamma, *a, *d);
            col1[j] = I1_inline(buf, t_full[k], t_full[j], gamma, *a, *d);
            col2[j] = I2_inline(buf, t_full[k], t_full[j], gamma, *a, *d);
        }
    } else {
        col0 = I0_loader->get_column(k);
        col1 = I1_loader->get_column(k);
        col2 = I2_loader->get_column(k);
    }

    arma::vec integrand = h0 % col0 + h1 % col1 + h2 % col2;

    double ic_val;
    if (use_inline_ic_integral) {
        ExponentialIC* eic = dynamic_cast<ExponentialIC*>(ic);
        if (eic) {
            ic_val = ic_int_inline(buf, t_full[k], gamma, *a, *d, x1 - buf, eic->lambda);
        } else {
            ic_val = 0.0;
        }
    } else {
        ic_val = IC_integral(k);
    }

    return g(k)*std::exp(-buf*buf) + trapezoidal_integrate(t_full, integrand, 0, t_curr, dt) + ic_val;
}

// -----------------------------------------------------------------------
// Main simulation loop (unchanged)
// -----------------------------------------------------------------------

FisherBasicSolver::Result FisherBasicSolver::run(std::ostream &logger, bool verbose) {
    arma::mat u(Nx, Nt);
    u.col(0) = ic->func(x);

    arma::vec u_curr = u.col(0);
    double t_curr = t0;

    init_matrix();
    b = arma::zeros<arma::vec>(Nx);

    if (verbose) {
        logger << "Beginning implicit finite differences simulation." << std::endl;
        logger << "Number of iterations: " << Nt*tout << "." << std::endl;
    }
    size_t output_num_iter = Nt / 20;

    logger << "u size: " << u.n_rows << " x " << u.n_cols << std::endl;
    logger << "t size: " << t.n_elem << std::endl;
    logger << "x size: " << x.n_elem << std::endl;
    logger << "u_curr size: " << u_curr.n_elem << std::endl;
    logger << "Matrix M size: " << M.n_rows << " x " << M.n_cols << std::endl;
    logger << "Vector b size: " << b.n_rows << " x " << std::endl;

    arma::vec g  = arma::zeros(t_full.n_elem);
    arma::vec gp = arma::zeros(t_full.n_elem);
    size_t idx_extract_g = (Nx-1) - (size_t)round(buf/dx);
    g(0) = u(idx_extract_g, 0);
    if (use_td_bc) {
        logger << "Initial g value: " << g(0) << std::endl;
    }

    size_t k = 0;
    for (size_t i = 0; i < Nt - 1; i += 1) {
        for (size_t j = 0; j < tout; j += 1) {
            update_rhs(u_curr, t_curr);

            if (use_td_bc) {
                if (k > 0) {
                    g(k) = u_curr(idx_extract_g);
                    gp(k) = (g(k) - g(k-1))/dt;
                }
                b(Nx-1) = compute_boundary_value(t_curr, k, g, gp);
            } else {
                b(Nx-1) = 0;
            }

            t_curr += dt;
            k += 1;

            update_matrix(t_curr);
            u_curr = arma::spsolve(M, b);
        }
        u.col(i+1) = u_curr;

        if (verbose && (i+1) % output_num_iter == 0) {
            logger << "Completed iteration " << (i+1)*tout << " (" << ((i+1)*100/Nt) << "%)." << std::endl;
            if (use_td_bc) {
                logger << "Current g value: " << g(k-1) << std::endl;
                logger << "Current gp value: " << gp(k-1) << std::endl;
                logger << "Current boundary value: " << b(Nx-1) << std::endl;
            }
        }
    }

    if (verbose) {
        logger << "Simulation complete." << std::endl;
    }

    Result sol;
    sol.t = t;
    sol.x = x;
    sol.u = u;

    return sol;
}
