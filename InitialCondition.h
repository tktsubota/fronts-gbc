#include <armadillo>

#pragma once

class InitialCondition {
public:
    virtual arma::vec func(const arma::vec &x) const = 0;
    
    virtual ~InitialCondition() {}
};