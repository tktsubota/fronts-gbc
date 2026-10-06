#include "options.h"

arma::vec GaussianIC::func(const arma::vec &x) const {
    return (1.0 / (sqrt(2.0 * M_PI) * sigma)) * arma::exp(-arma::square(x) / (2.0 * sigma*sigma));
}

arma::vec ExponentialIC::func(const arma::vec &x) const {
    return arma::exp(-lambda*arma::abs(x));
}

arma::vec TanhIC::func(const arma::vec &x) const {
    return 0.5*arma::tanh(-a*x) + 0.5;
}

arma::vec StepIC::func(const arma::vec& x) const {
    // 1.0 where x < 0, else 0.0
    return arma::conv_to<arma::vec>::from(x < 0);
}

arma::vec LineIC::func(const arma::vec& x) const {
    return arma::ones<arma::vec>(x.n_elem);
}


double ConstantCoefficient::func(double t) const {
    return c;
}

double LinearCoefficient::func(double t) const {
    return a + b*t;
}

double QuadraticCoefficient::func(double t) const {
    return a + b*t + c*t*t;
}

double ExponentialCoefficient::func(double t) const {
    return a*std::exp(b*t);
}

double AlgebraicCoefficient::func(double t) const {
    return a + b*std::pow(t,p);
}

double TanhCoefficient::func(double t) const {
    return a*std::tanh(t - b) + c;
}

double ConstantCoefficient::integral(double t) const {
    return c * t;
}

double LinearCoefficient::integral(double t) const {
    return a*t + (b/2.0)*t*t;
}

double QuadraticCoefficient::integral(double t) const {
    return a*t + (b/2.0)*t*t + (c/3.0)*t*t*t;
}

double ExponentialCoefficient::integral(double t) const {
    return (a/b)*std::exp(b*t);
}

double AlgebraicCoefficient::integral(double t) const {
    if (p == -1.0) {
        return a*t + b*std::log(t);
    } else {
        return a*t + b/(p+1.0)*std::pow(t, p+1.0);
    }
}

double TanhCoefficient::integral(double t) const {
    return a*std::log(std::cosh(t - b)) + c*t;
}