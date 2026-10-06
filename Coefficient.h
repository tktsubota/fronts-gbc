#include <armadillo>

#pragma once

class Coefficient {
public:
    virtual double func(double t) const = 0;
    virtual double integral(double t) const = 0; // indefinite integral; use integral(t)-integral(s) for ∫_s^t

    virtual ~Coefficient() {}
};