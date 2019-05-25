#include <vector> /* representing state */
#include <algorithm> /* for numeric min/max */
#include <cmath> /* exp */
#include <iostream> /* cout */

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

typedef double icing_float;
// Implementation borrowed from here: https://www.boost.org/doc/libs/1_53_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/tutorial/harmonic_oscillator.html
typedef std::vector<icing_float> state_type;

class ChaseIcing {
    // Each enum maps to an index within the state variable for accessing that
    // function's value.
    enum Function { BG_t = 0
        , Q_t
        , I_t
        , P1_t
        , P2_t
        , P_t
        , u_en_t
    };

    /* a bunch of constants that should be in the namespace */
    const icing_float n_I = 0.0;
    const icing_float d1 = 0.0;
    const icing_float d2 = 0.0;
    const icing_float P_max = 0.0;

    icing_float P_min(const auto &x, auto t)
    {
        return std::min<icing_float>(d2 * x[P2_t], P_max);
    }

    icing_float alpha_decay(const icing_float variable, const icing_float decay_parameter)
    {
        return variable / (1.0 + decay_parameter * variable);
    }

    icing_float Q_frac(const state_type &x, const icing_float t)
    {
        const icing_float a_G = 0.0;
        const auto Q = x[Q_t];
        return alpha_decay(Q, a_G);
    }

    // PN(t) -> Parenteral nutrition input, eg IV dextrose
    icing_float _P(const state_type &x, const icing_float t)
    {
        const auto _PN_t = 0.0; // TODO -> derive this over the network

        return P_min(x, t) + _PN_t;
    }

    icing_float BG_dot(const state_type &x, const icing_float t)
    {
        const icing_float p_G = 0.0;
        const icing_float S_I = 0.0;
        const icing_float BG = x[BG_t];
        const icing_float Q = x[Q_t];
        const icing_float P = _P(x, t);
        const icing_float EGP_b = 0.0;
        const icing_float CNS = 0.0;
        const icing_float V_G = 1.0;

        icing_float dBGdt = 0.0;
        // BG' = -p_G BG(t)
        //	 - S_I BG(t) (Q(t) / (1 + a_G Q(t)))
        //	 + (P(t) + EGP_b -CNS)/V_G
        dBGdt += -p_G * BG;
        dBGdt += -S_I*BG * Q_frac(x, t);
        dBGdt += (P + EGP_b - CNS)/V_G;

        return dBGdt;
    }

    icing_float Q_dot(const state_type &x, const icing_float t)
    {
        const icing_float I = x[I_t];
        const icing_float Q = x[Q_t];
        const icing_float n_C = 0.0;

        return n_I * (I - Q) - n_C * Q_frac(x, t);
    }

    auto _u_en(const auto &x, auto t)
    {
        const icing_float k1 = 0.0;
        const icing_float k2 = 0.0;
        const icing_float k3 = 1.0;
        return k1 * std::exp(-x[I_t] * (k2 / k3));
    }

    auto I_dot(const auto &x, auto t)
    {
        const icing_float n_K = 0.0;
        const icing_float I = x[I_t];
        const icing_float n_L = 0.0;
        const icing_float a_I = 0.0;
        const icing_float Q = x[Q_t];
        const icing_float u_ex = 0.0; // TODO: get this over the network
        const icing_float V_I = 1.0;
        const icing_float x_L = 0.0;
        const icing_float u_en = _u_en(x, t);

        auto dIdt = icing_float(0.0);
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
        const icing_float D = 0.0; // enteral feed rate TODO: get from network
        return (icing_float) -d1 * x[P1_t] + D;
    }

    auto P2_dot(const auto &x, auto t)
    {
        return -P_min(x, t) + d1 * x[P1_t];
    }

public:
    void operator() (const state_type &x, state_type &dxdt, const icing_float t)
    {
        dxdt[BG_t] = BG_dot(x, t);
        dxdt[Q_t] = Q_dot(x, t);
        dxdt[I_t] = I_dot(x, t);
        dxdt[P1_t] = P1_dot(x, t);
        dxdt[P2_t] = P2_dot(x, t);

        std::cout << "BG: " << dxdt[BG_t] << std::endl;
    }
};

int main(void)
{
    state_type x = {0.0, 0.0, 0.0, 0.0, 0.0};
    ChaseIcing model;
    const double dt = 0.1;
    runge_kutta4<state_type> stepper;
    integrate_const<runge_kutta4<state_type>, ChaseIcing>(stepper, model, x, 0.0, 10.0, dt);
	return 0;
}
