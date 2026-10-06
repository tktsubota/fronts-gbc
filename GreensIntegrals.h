#pragma once

#include <cmath>
#include "Coefficient.h"

/*
 * Closed-form Green's function integrals for the GBC method.
 *
 * I0, I1, I2, and ic_int are evaluated from their closed-form expressions.
 * The expression structure is chosen for numerical stability — do not merge
 * or simplify expressions.
 *
 * All functions take:
 *   x       - buffer position (buf)
 *   t_val   - current time
 *   s       - source time
 *   gamma   - frame speed multiplier
 *   a_coeff - Coefficient for a(t)
 *   d_coeff - Coefficient for d(t)
 *
 * Di(t,s) = d_coeff.integral(t) - d_coeff.integral(s)  = ∫_s^t d(t') dt'
 * Ai(t,s) = a_coeff.integral(t) - a_coeff.integral(s)  = ∫_s^t a(t') dt'
 *
 * All return 0 when t_val <= s (causality).
 */

inline double I0_inline(double x, double t_val, double s,
                        double gamma,
                        const Coefficient& a_coeff,
                        const Coefficient& d_coeff)
{
    if (t_val <= s) return 0.0;

    double Di_ts     = d_coeff.integral(t_val) - d_coeff.integral(s);
    double Ai_ts     = a_coeff.integral(t_val) - a_coeff.integral(s);
    double sqrt_term = std::sqrt(4.0 + 1.0 / Di_ts);
    double arg1      = (-x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);
    double arg2      = ( x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);
    double argouter  = -((x * (gamma + x)
                          - Ai_ts * (1.0 + 4.0 * Di_ts)
                          + gamma * Di_ts * (2.0 * x + gamma * Di_ts))
                         / (1.0 + 4.0 * Di_ts));
    double arginner  = (gamma * x) / (1.0 + 4.0 * Di_ts);

    double exp_outer = std::exp(argouter);

    double bracket = exp_outer * (-1.0 - std::erf(arg1))
                   + std::exp(argouter + arginner) * (1.0 + std::erf(arg2));

    return bracket / (2.0 * sqrt_term * std::sqrt(Di_ts));
}

inline double I1_inline(double x, double t_val, double s,
                        double gamma,
                        const Coefficient& a_coeff,
                        const Coefficient& d_coeff)
{
    if (t_val <= s) return 0.0;

    double Di_ts      = d_coeff.integral(t_val) - d_coeff.integral(s);
    double Ai_ts      = a_coeff.integral(t_val) - a_coeff.integral(s);
    double sqrt_term  = std::sqrt(4.0 + 1.0 / Di_ts);

    double outer_exp_term = Ai_ts - (x + gamma * Di_ts) * (x + gamma * Di_ts) / (4.0 * Di_ts);

    double expA = std::exp(outer_exp_term
                  + (gamma - x / Di_ts) * (gamma - x / Di_ts) / (4.0 * (4.0 + 1.0 / Di_ts)));
    double expB = std::exp(outer_exp_term
                  + (x + gamma * Di_ts) * (x + gamma * Di_ts) / (4.0 * Di_ts * (1.0 + 4.0 * Di_ts)));

    double arg1 = (-x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);
    double arg2 = ( x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);

    double bracket = (
        -gamma * Di_ts * (
            expA
            - expB
            + expA * std::erf(arg1)
            - expB * std::erf(arg2)
        )
        + x * (
            expA
            + expB
            + expA * std::erf(arg1)
            + expB * std::erf(arg2)
        )
    );

    double denom = 2.0 * sqrt_term * std::sqrt(Di_ts) * (1.0 + 4.0 * Di_ts);
    return bracket / denom;
}

inline double I2_inline(double x, double t_val, double s,
                        double gamma,
                        const Coefficient& a_coeff,
                        const Coefficient& d_coeff)
{
    if (t_val <= s) return 0.0;

    double Di_ts      = d_coeff.integral(t_val) - d_coeff.integral(s);
    double Ai_ts      = a_coeff.integral(t_val) - a_coeff.integral(s);
    double sqrt_term  = std::sqrt(4.0 + 1.0 / Di_ts);
    double sqrt_pi    = std::sqrt(M_PI);

    double outer_exp_term = Ai_ts - (x + gamma * Di_ts) * (x + gamma * Di_ts) / (4.0 * Di_ts);

    double expA = std::exp(outer_exp_term
                  + (gamma - x / Di_ts) * (gamma - x / Di_ts) / (4.0 * (4.0 + 1.0 / Di_ts)));
    double expB = std::exp(outer_exp_term
                  + (x + gamma * Di_ts) * (x + gamma * Di_ts) / (4.0 * Di_ts * (1.0 + 4.0 * Di_ts)));

    double arg1 = (-x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);
    double arg2 = ( x + gamma * Di_ts) / (2.0 * sqrt_term * Di_ts);

    double common_bracket = (
        expA
        - expB
        + expA * std::erf(arg1)
        - expB * std::erf(arg2)
    );

    double bracket = (
        -sqrt_pi * x * x * common_bracket
        - (8.0 + gamma * gamma) * sqrt_pi * Di_ts * Di_ts * common_bracket
        + 2.0 * Di_ts * (
            -expA * sqrt_pi
            + expB * sqrt_pi
            + expA * gamma * sqrt_pi * x
            + expB * gamma * sqrt_pi * x
            + std::exp(outer_exp_term) * 2.0 * x * sqrt_term
            + expA * sqrt_pi * (-1.0 + gamma * x) * std::erf(arg1)
            + expB * sqrt_pi * ( 1.0 + gamma * x) * std::erf(arg2)
        )
    );

    double denom = 2.0 * sqrt_pi * sqrt_term * std::sqrt(Di_ts) * (1.0 + 4.0 * Di_ts) * (1.0 + 4.0 * Di_ts);
    return bracket / denom;
}

/*
 * Integral of the Green's function against an exponential initial condition
 * u(x,0) = exp(-lda * |x|) evaluated at x=buf for time t_val.
 *
 * Arguments:
 *   x        = buf
 *   boundary = x1 - buf
 *   lda      = IC steepness
 *
 * Returns 0 when t_val <= 0.
 */
inline double ic_int_inline(double x, double t_val,
                            double gamma,
                            const Coefficient& a_coeff,
                            const Coefficient& d_coeff,
                            double boundary, double lda)
{
    if (t_val <= 0.0) return 0.0;

    double Di_t0 = d_coeff.integral(t_val) - d_coeff.integral(0.0);
    double Ai_t0 = a_coeff.integral(t_val) - a_coeff.integral(0.0);

    double sqrt_inv_D  = std::sqrt(1.0 / Di_t0);
    double sqrt_term_4 = std::sqrt(4.0 + 1.0 / Di_t0);
    double denom_4     = 1.0 + 4.0 * Di_t0;

    // term A
    double expA  = std::exp(-boundary * lda + Ai_t0
                            - (gamma - lda) * (x + lda * Di_t0));
    double argA  = 0.5 * sqrt_inv_D * (-x + (gamma - 2.0 * lda) * Di_t0);
    double termA = -expA * (1.0 + std::erf(argA)) / sqrt_inv_D;

    // term B
    double numB  = (x * (gamma + x - lda)
                    + (-lda * lda + gamma * (2.0 * x + lda)) * Di_t0
                    + gamma * gamma * Di_t0 * Di_t0);
    double expB  = std::exp(-boundary * lda + Ai_t0 - numB / denom_4);
    double argB  = (-x + (gamma - 2.0 * lda) * Di_t0) / (2.0 * sqrt_term_4 * Di_t0);
    double termB = expB * (1.0 + std::erf(argB)) / sqrt_term_4;

    // term C
    double expC  = std::exp(-boundary * lda + Ai_t0
                            - lda * (x + (gamma - lda) * Di_t0));
    double argC  = 0.5 * sqrt_inv_D * (x + (gamma - 2.0 * lda) * Di_t0);
    double termC = expC * (1.0 + std::erf(argC)) / sqrt_inv_D;

    // term D
    double numD  = (x * (x + lda)
                    + (-lda * lda + gamma * (2.0 * x + lda)) * Di_t0
                    + gamma * gamma * Di_t0 * Di_t0);
    double expD  = std::exp(-boundary * lda + Ai_t0 - numD / denom_4);
    double argD  = (x + (gamma - 2.0 * lda) * Di_t0) / (2.0 * sqrt_term_4 * Di_t0);
    double termD = expD * (-1.0 - std::erf(argD)) / sqrt_term_4;

    double bracket = termA + termB + termC + termD;
    return bracket / (2.0 * std::sqrt(Di_t0));
}
