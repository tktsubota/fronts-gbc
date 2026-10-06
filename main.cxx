#include <cmath>
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <H5Cpp.h>
#include <memory>
#include <chrono>
#include <map>
#include <stdexcept>

#include <tclap/CmdLine.h>
#include <string>
#include <vector>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>

#include "InitialCondition.h"
#include "Coefficient.h"
#include "FisherBasicSolver.h"
#include "Fisher23Solver.h"
#include "options.h"

std::string get_iso8601_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto now_time = std::chrono::system_clock::to_time_t(now);
    auto now_tm = *std::gmtime(&now_time); // Convert to UTC time

    std::ostringstream oss;
    oss << std::put_time(&now_tm, "%Y%m%dT%H%M%SZ");
    return oss.str();
}

std::vector<double> scale_params(const std::string& coeffType,
                                 const std::vector<double>& originalParams,
                                 double scaleFactor) {
    std::vector<double> scaled;
    if (coeffType == "const" || coeffType == "linear" || coeffType == "quadratic") {
        for (double p : originalParams) {
            scaled.push_back(scaleFactor * p);
        }
    } else if (coeffType == "exp") {
        if (originalParams.size() == 2) {
            scaled.push_back(scaleFactor * originalParams[0]);
            scaled.push_back(originalParams[1]);
        } else if (originalParams.size() == 1) {
            scaled.push_back(scaleFactor);
            scaled.push_back(originalParams[0]);
        } else {
            scaled.push_back(scaleFactor);
            scaled.push_back(1.0);
        }
    } else if (coeffType == "algebraic") {
        if (originalParams.size() == 3) {
            scaled.push_back(scaleFactor * originalParams[0]);
            scaled.push_back(scaleFactor * originalParams[1]);
            scaled.push_back(originalParams[2]);
        } else if (originalParams.size() == 2) {
            scaled.push_back(scaleFactor);
            scaled.push_back(scaleFactor * originalParams[0]);
            scaled.push_back(originalParams[1]);
        } else if (originalParams.size() == 1) {
            scaled.push_back(scaleFactor);
            scaled.push_back(scaleFactor);
            scaled.push_back(originalParams[0]);
        } else {
            scaled.push_back(scaleFactor);
            scaled.push_back(scaleFactor);
            scaled.push_back(0.5);
        }
    } else if (coeffType == "tanh") {
        if (originalParams.size() == 3) {
            scaled.push_back(scaleFactor * originalParams[0]);
            scaled.push_back(originalParams[1]);
            scaled.push_back(scaleFactor * originalParams[2]);
        } else if (originalParams.size() == 2) {
            scaled.push_back(scaleFactor * originalParams[0]);
            scaled.push_back(originalParams[1]);
        } else if (originalParams.size() == 1) {
            scaled.push_back(scaleFactor * originalParams[0]);
        } else {
            scaled.push_back(scaleFactor);
            scaled.push_back(1.0);
        }
    } else {
        throw std::runtime_error("scale_params: coeff type not implemented");
    }
    return scaled;
}

std::unique_ptr<Coefficient> make_coefficient(const std::string& coeffType,
                                             const std::vector<double>& params) {
    if (coeffType == "const") {
        if (params.size() == 1) return std::make_unique<ConstantCoefficient>(params[0]);
        return std::make_unique<ConstantCoefficient>(1.0);

    } else if (coeffType == "linear") {
        if (params.size() == 2) return std::make_unique<LinearCoefficient>(params[0], params[1]);
        if (params.size() == 1) return std::make_unique<LinearCoefficient>(1.0, params[0]);
        return std::make_unique<LinearCoefficient>(1.0, 1.0);

    } else if (coeffType == "quadratic") {
        if (params.size() == 3) return std::make_unique<QuadraticCoefficient>(params[0], params[1], params[2]);
        if (params.size() == 2) return std::make_unique<QuadraticCoefficient>(params[0], 0.0, params[1]);
        if (params.size() == 1) return std::make_unique<QuadraticCoefficient>(1.0, 0.0, params[0]);
        return std::make_unique<QuadraticCoefficient>(1.0, 0.0, 1.0);

    } else if (coeffType == "exp") {
        if (params.size() == 2) return std::make_unique<ExponentialCoefficient>(params[0], params[1]);
        if (params.size() == 1) return std::make_unique<ExponentialCoefficient>(1.0, params[0]);
        return std::make_unique<ExponentialCoefficient>(1.0, 1.0);

    } else if (coeffType == "algebraic") {
        if (params.size() == 3) return std::make_unique<AlgebraicCoefficient>(params[0], params[1], params[2]);
        if (params.size() == 2) return std::make_unique<AlgebraicCoefficient>(1.0, params[0], params[1]);
        if (params.size() == 1) return std::make_unique<AlgebraicCoefficient>(1.0, 1.0, params[0]);
        return std::make_unique<AlgebraicCoefficient>(1.0, 1.0, 0.5);

    } else if (coeffType == "tanh"){
        if (params.size() == 3) return std::make_unique<TanhCoefficient>(params[0], params[1], params[2]);
        if (params.size() == 2) return std::make_unique<TanhCoefficient>(params[0], params[1], 0.0);
        if (params.size() == 1) return std::make_unique<TanhCoefficient>(params[0], 0.0, 0.0);
        return std::make_unique<TanhCoefficient>(1.0, 0.0, 0.0); 

    } else {
        throw std::runtime_error("Unknown coefficient type: " + coeffType);
    }
}

struct ParsedCoefficients {
    std::map<std::string, std::unique_ptr<Coefficient>> coeffs;
    double gamma;
};

ParsedCoefficients parse_coefficients(const std::vector<std::string>& tokens) {
    std::map<std::string, std::string> coeffTypes;
    std::map<std::string, std::vector<double>> coeffParams;

    size_t i = 0;
    while (i < tokens.size()) {
        std::string coeffName = tokens[i];
        std::string coeffType = tokens[i + 1];
        std::vector<double> params;

        size_t j = i + 2;
        while (j < tokens.size() && tokens[j] != "a" && tokens[j] != "d" && tokens[j] != "gamma") {
            try {
                params.push_back(std::stod(tokens[j]));
            } catch (...) {}
            ++j;
        }

        coeffTypes[coeffName] = coeffType;
        coeffParams[coeffName] = params;
        i = j;
    }

    // defaults
    if (coeffTypes.find("a") == coeffTypes.end()) {
        coeffTypes["a"] = "const";
    }
    if (coeffTypes.find("d") == coeffTypes.end()) {
        coeffTypes["d"] = "const";
    }
    if (coeffTypes.find("gamma") == coeffTypes.end()) {
        coeffTypes["gamma"] = "manual";
    }

    // build a, d
    auto a = make_coefficient(coeffTypes["a"], coeffParams["a"]);
    auto d = make_coefficient(coeffTypes["d"], coeffParams["d"]);

    // gamma handling
    double gamma;
    if (coeffParams["gamma"].empty()) {
        gamma = 2.0; // default
    } else {
        gamma = coeffParams["gamma"][0];
    }

    // scale d to get c
    auto scaledParams = scale_params(coeffTypes["d"], coeffParams["d"], gamma);
    auto c = make_coefficient(coeffTypes["d"], scaledParams);

    ParsedCoefficients result;
    result.gamma = gamma;
    result.coeffs["a"] = std::move(a);
    result.coeffs["d"] = std::move(d);
    result.coeffs["c"] = std::move(c);

    return result;
}

int main(int argc, char* argv[]) {
    #ifdef ARMA_USE_SUPERLU
        std::cout << "SuperLU enabled!" << std::endl;
    #else
        std::cout << "SuperLU NOT enabled." << std::endl;
    #endif
    try {
        // Create the command-line parser
        TCLAP::CmdLine cmd("Fisher-KPP front simulator with Green's function boundary conditions", ' ', "1.0.0");

        // Option arguments for the PDE and initial condition
        TCLAP::ValueArg<std::string> pde_arg("p", "pde", "PDE name (e.g. 'fisher-basic', 'fisher23')", false, "fisher-basic", "string");
        TCLAP::ValueArg<std::string> pde_params_arg("P", "pdeParams", "PDE parameters", false, "a const 1 d const 1 gamma manual 2", "string");
        TCLAP::ValueArg<std::string> ic_arg("i", "ic", "Initial condition (e.g. 'gaussian')", true, "gaussian 0.5", "string");

        // Numerical arguments
        TCLAP::ValueArg<double> t0_arg("t", "t0", "Initial time (undefined behavior if > 0)", false, 0.0, "double");
        TCLAP::ValueArg<double> t1_arg("T", "t1", "Final time", false, 10.0, "double");
        TCLAP::ValueArg<double> dt_arg("d", "dt", "Time step", false, 0.1, "double");
        TCLAP::ValueArg<double> x0_arg("x", "x0", "Left spatial boundary", false, -20.0, "double");
        TCLAP::ValueArg<double> x1_arg("X", "x1", "Right spatial boundary", false, 20.0, "double");
        TCLAP::ValueArg<double> dx_arg("D", "dx", "Spatial step", false, 0.01, "double");
        TCLAP::ValueArg<int> tout_arg("o", "tout", "Output increment (# time steps before output)", false, 1, "int");
        TCLAP::ValueArg<double> buf_arg("b", "buf", "Time-dependent boundary condition buffer size (if zero, use naive Dirichlet; if -1, use naive Neumann)", false, 0, "double");
        TCLAP::ValueArg<std::string> impulse_arg("m", "impulse", "Impulse HDF5 file", false, "", "string");
        TCLAP::ValueArg<std::string> ic_integral_arg("n", "icintegral", "IC integral HDF5 file", false, "", "string");
        TCLAP::SwitchArg quiet_arg("q", "quiet", "Quiet output", false);

        // Add all arguments to the command-line parser
        cmd.add(pde_arg);
        cmd.add(pde_params_arg);
        cmd.add(ic_arg);
        cmd.add(t0_arg);
        cmd.add(t1_arg);
        cmd.add(dt_arg);
        cmd.add(x0_arg);
        cmd.add(x1_arg);
        cmd.add(dx_arg);
        cmd.add(tout_arg);
        cmd.add(buf_arg);
        cmd.add(impulse_arg);
        cmd.add(ic_integral_arg);
        cmd.add(quiet_arg);

        // Parse the arguments
        cmd.parse(argc, argv);

        // Helper to split space-separated string
        auto split = [](const std::string &s) -> std::vector<std::string> {
            std::vector<std::string> tokens;
            std::istringstream iss(s);
            std::string token;
            while (iss >> token) {
                tokens.push_back(token);
            }
            return tokens;
        };

        // Access the parsed arguments
        std::string pde_str = pde_arg.getValue();
        std::string pde_params_str = pde_params_arg.getValue();
        std::vector<std::string> pde_params = split(pde_params_str);

        std::string ic_str = ic_arg.getValue();
        std::vector<std::string> ic_params = split(ic_str);

        double t0 = t0_arg.getValue();
        double t1 = t1_arg.getValue();
        double dt = dt_arg.getValue();
        double x0 = x0_arg.getValue();
        double x1 = x1_arg.getValue();
        double dx = dx_arg.getValue();
        int tout = tout_arg.getValue();
        double buf = buf_arg.getValue();
        std::string impulse_file = impulse_arg.getValue();
        std::string ic_integral_file = ic_integral_arg.getValue();
        bool quiet = quiet_arg.getValue();

        // Output the parsed arguments
        std::cout << "---" << std::endl;
        std::cout << "PDE: " << pde_str << " (";
        for (const auto &param: pde_params) {
            std::cout << " " << param;
        }
        std::cout << " )" << std::endl;

        std::cout << "Initial Condition: " << ic_str << std::endl;

        std::cout << "t0: " << t0 << ", t1: " << t1 << ", dt: " << dt << std::endl;
        std::cout << "x0: " << x0 << ", x1: " << x1 << ", dx: " << dx << std::endl;
        std::cout << "tout: " << tout << std::endl;
        std::cout << "buf: " << buf << std::endl;
        std::cout << "---" << std::endl;

        // Create output directory if it does not exist
        if (mkdir("sims", 0755) != 0 && errno != EEXIST) {
            std::cerr << "Error: Could not create directory sims: " << std::strerror(errno) << std::endl;
            return 1;
        }

        // Create the coefficients
        ParsedCoefficients parsed_coeffs = parse_coefficients(pde_params);

        // Create the initial condition
        std::unique_ptr<InitialCondition> ic;
        if (ic_params[0] == "gaussian") {
            double sigma = (ic_params.size() > 1) ? std::stod(ic_params[1]) : 1.0;
            ic = std::make_unique<GaussianIC>(sigma);
        } else if ((ic_params[0] == "exp") || (ic_params[0] == "exponential")) {
            double lambda = (ic_params.size() > 1) ? std::stod(ic_params[1]) : 1.0;
            ic = std::make_unique<ExponentialIC>(lambda);
        } else if (ic_params[0] == "tanh") {
            double a = (ic_params.size() > 1) ? std::stod(ic_params[1]) : 1.0;
            ic = std::make_unique<TanhIC>(a);
        } else if (ic_params[0] == "step") {
            ic = std::make_unique<StepIC>();
        } else if (ic_params[0] == "line") {
            ic = std::make_unique<LineIC>();
        }

        // Create the solver
        double gamma = parsed_coeffs.gamma;
        std::unique_ptr<FisherBasicSolver> solver;
        if (pde_str == "fisher23") {
            std::cout << "Using Fisher23 Solver." << std::endl;
            solver = std::make_unique<Fisher23Solver>(parsed_coeffs.coeffs["a"].get(), parsed_coeffs.coeffs["d"].get(), parsed_coeffs.coeffs["c"].get(), t0, t1, dt, x0, x1, dx, ic.get(), tout, buf, gamma, impulse_file, ic_integral_file);
        } else {
            std::cout << "Using Fisher Basic Solver." << std::endl;
            solver = std::make_unique<FisherBasicSolver>(parsed_coeffs.coeffs["a"].get(), parsed_coeffs.coeffs["d"].get(), parsed_coeffs.coeffs["c"].get(), t0, t1, dt, x0, x1, dx, ic.get(), tout, buf, gamma, impulse_file, ic_integral_file);
        }
        // Run the simulation
        auto start = std::chrono::high_resolution_clock::now();

        FisherBasicSolver::Result test = solver->run(std::cout, !quiet);

        // Compute time
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_sec = std::chrono::duration<double>(end - start).count();
        int hours = static_cast<int>(elapsed_sec) / 3600;
        int minutes = (static_cast<int>(elapsed_sec) % 3600) / 60;
        double seconds = elapsed_sec - hours * 3600 - minutes * 60;
        std::cout << "Simulation time: ";
        if (hours > 0) std::cout << hours << " h ";
        if (minutes > 0 || hours > 0) std::cout << minutes << " m ";
        std::cout << seconds << " s" << std::endl;

        // Save the simulation result
        // Generate filename with timestamp
        std::string end_timestamp = get_iso8601_timestamp();
        std::string filename = "sims/" + end_timestamp + "_sim" + ".h5";

        H5::H5File file(filename, H5F_ACC_TRUNC);

        // Save time vector
        H5::DataSpace t_dataspace(1, (hsize_t[]){(hsize_t)test.t.n_elem});
        H5::DataSet t_dataset = file.createDataSet("t", H5::PredType::NATIVE_DOUBLE, t_dataspace);
        t_dataset.write(test.t.memptr(), H5::PredType::NATIVE_DOUBLE);

        // Save x vector
        H5::DataSpace x_dataspace(1, (hsize_t[]){(hsize_t)test.x.n_elem});
        H5::DataSet x_dataset = file.createDataSet("x", H5::PredType::NATIVE_DOUBLE, x_dataspace);
        x_dataset.write(test.x.memptr(), H5::PredType::NATIVE_DOUBLE);

        // Save u matrix
        H5::DataSpace u_dataspace(2, (hsize_t[]){(hsize_t)test.u.n_cols, (hsize_t)test.u.n_rows});
        H5::DataSet u_dataset = file.createDataSet("u", H5::PredType::NATIVE_DOUBLE, u_dataspace);
        u_dataset.write(test.u.memptr(), H5::PredType::NATIVE_DOUBLE);

        // Variable-length string type
        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);

        // Store PDE string and parameters
        file.createAttribute("pde", H5::StrType(H5::PredType::C_S1, pde_str.size()), H5::DataSpace()).write(H5::StrType(H5::PredType::C_S1, pde_str.size()), pde_str.c_str());
        file.createAttribute("pde_params", H5::StrType(H5::PredType::C_S1, pde_params_str.size()), H5::DataSpace()).write(H5::StrType(H5::PredType::C_S1, pde_params_str.size()), pde_params_str.c_str());
        file.createAttribute("ic", H5::StrType(H5::PredType::C_S1, ic_str.size()), H5::DataSpace()).write(H5::StrType(H5::PredType::C_S1, ic_str.size()), ic_str.c_str());

        // Store metadata attributes
        file.createAttribute("t0", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE,  &t0);
        file.createAttribute("t1", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &t1);
        file.createAttribute("dt", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &dt);
        file.createAttribute("x0", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &x0);
        file.createAttribute("x1", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &x1);
        file.createAttribute("dx", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &dx);
        file.createAttribute("tout", H5::PredType::NATIVE_INT, H5::DataSpace()).write(H5::PredType::NATIVE_INT, &tout);
        file.createAttribute("buf", H5::PredType::NATIVE_DOUBLE, H5::DataSpace()).write(H5::PredType::NATIVE_DOUBLE, &buf);

        // Save simulation file
        file.close();
        std::cout << "Saved simulation to file: " << filename << std::endl;

        // Save the simulation parameters to csv
        std::string sims_csv_filename = "sims/sims.csv";
        bool sims_csv_exists = std::ifstream(sims_csv_filename).good();
        std::ofstream sims_csv(sims_csv_filename, std::ios::app);
        if (!sims_csv.is_open()) {
            std::cerr << "Error: Could not open file " << sims_csv_filename << std::endl;
            return 1;
        }

        // Write header if the CSV is new
        if (!sims_csv_exists) {
            sims_csv << "timestamp,pde,pde_params,ic,t0,t1,dt,x0,x1,dx,buf,tout\n";
        }

        // Append data to CSV
        sims_csv << end_timestamp << ","
            << "\"" << pde_str << "\","
            << "\"" << pde_params_str << "\","
            << "\"" << ic_str << "\","
            << t0 << "," << t1 << "," << dt << ","
            << x0 << "," << x1 << "," << dx << ","
            << buf << ","
            << tout << "\n";

        sims_csv.close();

    } catch (TCLAP::ArgException &e) {
        // Catch any argument parsing errors
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        return 1;
    }

    return 0;
}