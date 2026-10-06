#include <cmath>
#include <armadillo>
#include "InitialCondition.h"
#include "Coefficient.h"

class GaussianIC: public InitialCondition {
private:
    double sigma = 1.0;
public:
    GaussianIC(double sigma): sigma(sigma) {};
    arma::vec func(const arma::vec &x) const override;
};

class ExponentialIC: public InitialCondition {
public:
    double lambda = 1.0;
    ExponentialIC(double lambda): lambda(lambda) {};
    arma::vec func(const arma::vec &x) const override;
};

class TanhIC: public InitialCondition {
private:
    double a = 1.0;
public:
    TanhIC(double a): a(a) {};
    arma::vec func(const arma::vec &x) const override;
};


class StepIC: public InitialCondition {
public:
    arma::vec func(const arma::vec &x) const override;
};


class LineIC: public InitialCondition {
public:
    arma::vec func(const arma::vec &x) const override;
};


class ConstantCoefficient: public Coefficient {
private:
    double c = 1.0;
public:
    ConstantCoefficient(double c): c(c) {};
    double func(double t) const override;
    double integral(double t) const override;
};

class LinearCoefficient: public Coefficient {
private:
    double a = 1.0;
    double b = 1.0;
public:
    LinearCoefficient(double a, double b): a(a), b(b) {};
    double func(double t) const override;
    double integral(double t) const override;
};

class QuadraticCoefficient: public Coefficient {
private:
    double a = 1.0;
    double b = 0.0;
    double c = 1.0;
public:
    QuadraticCoefficient(double a, double b, double c): a(a), b(b), c(c) {};
    double func(double t) const override;
    double integral(double t) const override;
};

class ExponentialCoefficient: public Coefficient {
private:
    double a = 1.0;
    double b = 1.0;
public:
    ExponentialCoefficient(double a, double b): a(a), b(b) {};
    double func(double t) const override;
    double integral(double t) const override;
};


class AlgebraicCoefficient: public Coefficient {
private:
    double a = 1.0;
    double b = 1.0;
    double p = 0.5;
public:
    AlgebraicCoefficient(double a, double b, double p): a(a), b(b), p(p) {};
    double func(double t) const override;
    double integral(double t) const override;
};

class TanhCoefficient: public Coefficient {
private:
    double a = 1.0;
    double b = 0.0;
    double c = 0.0;
public:
    TanhCoefficient(double a, double b, double c): a(a), b(b), c(c) {};
    double func(double t) const override;
    double integral(double t) const override;
};