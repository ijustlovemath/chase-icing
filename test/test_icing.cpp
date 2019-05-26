#include <vector> /* representing state */
#include <algorithm> /* for numeric min/max */
#include <cmath> /* exp, fmod */
#include <iostream> /* cout */

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

typedef double icing_float;
typedef std::vector<icing_float> state_type;

// [1] - Reference Paper: TBME-01160-2016 
// [2] - Reference Paper: computer methods and programs in biomedicine 102 (2011) 192–205
template<typename U = icing_float>
class ChaseIcing {
    // Each enum maps to an index within the state variable for accessing that
    // function's value.
    enum Function { G_t = 0
        , Q_t
        , I_t
        , P1_t
        , P2_t
        , P_t
        , u_en_t
    };

    /* a bunch of constants that should be in the namespace */
    const U n_I = 0.0075; // 1/min, Table II, [1]
    const U d1 = -std::log(0.5)/20.0; // 1/min, Table II, [1]
    const U d2 = -std::log(0.5)/100.0; // 1/min, Table II, [1]
    const U P_max = 6.11; // mmol/min, Table II, [1]

    auto P_min(const auto &x, auto t)
    {
        return std::min(d2 * x[P2_t], P_max);
    }

    auto alpha_decay(const auto variable, const auto decay_parameter)
    {
        return variable / (1.0 + decay_parameter * variable);
    }

    auto Q_frac(const auto &x, const auto t)
    {
        const U a_G = 1.0 / 65.0; // Table II, [1]
        const U Q = x[Q_t];
        return alpha_decay(Q, a_G);
    }

    // PN(t) -> Parenteral nutrition input, eg IV dextrose
    auto _P(const auto &x, const auto t)
    {
        const U _PN_t = 0.0; // TODO -> derive this over the network

        return P_min(x, t) + _PN_t;
    }

    auto _u_en(const auto &x, const auto t)
    {
        const U k1 = 14.9; // mU * l / mmol/min, Table II, [1]
        const U k2 = -49.9; // mU/min, Table II, [1]
        const U u_min = 16.7; // mU/min, Table II, [1]
        const U u_max = 266.7; // mU/min, Table II, [1]
        const U G = x[G_t];
        return std::min(std::max(u_min, k1 * G + k2), u_max);
    }

    auto G_dot(const auto &x, const auto t)
    {
        const U p_G = 0.006; // End of section 4.1, in [2]
        const U S_I = 0.5e-3; // TODO: patient specific
        const U G = x[G_t];
        const U Q = x[Q_t];
        const U P = _P(x, t);
        const U EGP = 1.16; // mmol/min Table II, [1]
        const U CNS = 0.3; // mmol/min Table II, [1]
        const U V_G = 13.3; // L, Table II, [1]

        U dGdt = 0.0;
        // G' = -p_G G(t)
        //	 - S_I G(t) (Q(t) / (1 + a_G Q(t)))
        //	 + (P(t) + EGP -CNS)/V_G
        dGdt += -p_G * G;
        dGdt += -S_I*G * Q_frac(x, t);
        dGdt += (P + EGP - CNS)/V_G;

        return dGdt;
    }

    auto Q_dot(const auto &x, const auto t)
    {
        const U I = x[I_t];
        const U Q = x[Q_t];
        const U n_C = 0.0075; // 1/min, Table II, [1]

        return n_I * (I - Q) - n_C * Q_frac(x, t);
    }

    auto I_dot(const auto &x, auto t)
    {
        const U n_K = 0.0542; // 1/min, Table II, [1]
        const U n_L = 0.1578; // 1/min, Table II, [1] 
        const U a_I = 1.7e-3; // 1/mU, Table II, [1]
        const U V_I = 4.0; // L, Table II, [1]
        const U x_L = 0.67; // unitless, Table II, [1]

        const U u_ex = 0.0; // TODO: get this over the network

        const U Q = x[Q_t];
        const U I = x[I_t];
        const U u_en = _u_en(x, t);

        auto dIdt = U(0.0);
        // I' = - n_K I(t)
        //	- n_L (I(t)/(1+a_I I(t)))
        //	- n_I (I(t) - Q(t))
        //	+ u_ex(t) / V_I
        //	+ (1 - x_L) u_en / V_I
        dIdt += -n_K * I;
        dIdt += -n_L * alpha_decay(I, a_I);
        dIdt += -n_I * (I - Q);
        dIdt += u_ex / V_I;
        dIdt += (1.0 - x_L) * u_en / V_I;
        return dIdt;
    }

    auto P1_dot(const auto &x, auto t)
    {
        const auto D = U(0.0); // enteral feed rate TODO: get from network
        return -d1 * x[P1_t] + D;
    }

    auto P2_dot(const auto &x, auto t)
    {
        return -P_min(x, t) + d1 * x[P1_t];
    }

public:
    using state_type = std::vector<U>;
    void operator() (const state_type &x, state_type &dxdt, const U t)
    {
        dxdt[G_t] = G_dot(x, t);
        dxdt[Q_t] = Q_dot(x, t);
        dxdt[I_t] = I_dot(x, t);
        dxdt[P1_t] = P1_dot(x, t);
        dxdt[P2_t] = P2_dot(x, t);

        if(std::fmod(t, 10.0) <= 0.1) {
            std::cout << 
            "t: " << t << "\t"
            "G: " << x[G_t] << std::endl;
        }
    }
};

int main(void)
{
    using Chase = ChaseIcing<>;
    using vec5 = Chase::state_type;

    const double dt = 0.1;

    vec5 x = {10.0, 3.0, 10.0, 4.0, 3.0};
    runge_kutta4<vec5> stepper;
    Chase model;  
    
    integrate_const(stepper, model, x, 0.0, 1000.0, dt);
	
    return 0;
}
