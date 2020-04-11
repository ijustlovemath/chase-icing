#include <vector> /* representing state */
#include <algorithm> /* for numeric min/max */
#include <cmath> /* log, exp, fmod */

#include <boost/numeric/odeint.hpp>
#include <ios>
using namespace boost::numeric::odeint;

typedef double icing_float;

// [1] - Reference Paper: TBME-01160-2016 
// [2] - Reference Paper: computer methods and programs in biomedicine 102 (2011) 192–205
template<typename U = icing_float>
class ChaseIcing {
    // Each enum maps to an index within the state vector for accessing that
    // function's value.
    enum Function { fn_G = 0
        , fn_Q
        , fn_I
        , fn_P1
        , fn_P2
        , fn_P
        , fn_u_en
    };
    using _state_type = std::vector<U>;

    /* a bunch of constants that should be in the namespace */
    const U n_I = 0.0075; // 1/min, Table II, [1]
    const U d1 = -std::log(0.5)/20.0; // 1/min, Table II, [1]
    const U d2 = -std::log(0.5)/100.0; // 1/min, Table II, [1]
    const U P_max = 6.11; // mmol/min, Table II, [1]

    auto P_min(const _state_type &x)
    {
        return std::min(d2 * x[fn_P2], P_max);
    }

    auto alpha_decay(const U variable, const U decay_parameter)
    {
        return variable / (1.0 + decay_parameter * variable);
    }

    auto Q_frac(const _state_type &x)
    {
        const U a_G = 1.0 / 65.0; // Table II, [1]
        const U Q = x[fn_Q];
        return alpha_decay(Q, a_G);
    }

    // PN(t) -> Parenteral nutrition input, eg IV dextrose
    auto _P(const _state_type &x)
    {
        const U _PN_ext = dextrose_rate; // TODO -> derive this over the network

        return P_min(x) + _PN_ext;
    }

    auto _u_en(const _state_type &x)
    {
        const U k1 = 14.9; // mU * l / mmol/min, Table II, [1]
        const U k2 = -49.9; // mU/min, Table II, [1]
        const U u_min = 16.7; // mU/min, Table II, [1]
        const U u_max = 266.7; // mU/min, Table II, [1]
        const U G = x[fn_G];
        return std::min(std::max(u_min, k1 * G + k2), u_max);
    }

    auto G_dot(const _state_type &x)
    {
        const U p_G = 0.006; // End of section 4.1, in [2]
        const U S_I = 0.5e-3; // TODO: patient specific
        const U G = x[fn_G];
        const U Q = x[fn_Q];
        const U P = _P(x);
        const U EGP = 1.16; // mmol/min Table II, [1]
        const U CNS = 0.3; // mmol/min Table II, [1]
        const U V_G = 13.3; // L, Table II, [1]

        U dGdt = 0.0;
        // G' = -p_G G(t)
        //	 - S_I G(t) (Q(t) / (1 + a_G Q(t)))
        //	 + (P(t) + EGP -CNS)/V_G
        dGdt += -p_G * G;
        dGdt += -S_I*G * Q_frac(x);
        dGdt += (P + EGP - CNS)/V_G;

        return dGdt;
    }

    auto Q_dot(const _state_type &x)
    {
        const U I = x[fn_I];
        const U Q = x[fn_Q];
        const U n_C = 0.0075; // 1/min, Table II, [1]

        return n_I * (I - Q) - n_C * Q_frac(x);
    }

    auto I_dot(const _state_type &x)
    {
        const U n_K = 0.0542; // 1/min, Table II, [1]
        const U n_L = 0.1578; // 1/min, Table II, [1] 
        const U a_I = 1.7e-3; // 1/mU, Table II, [1]
        const U V_I = 4.0; // L, Table II, [1]
        const U x_L = 0.67; // unitless, Table II, [1]

        const U u_ex = insulin_rate; // TODO: get this over the network

        const U Q = x[fn_Q];
        const U I = x[fn_I];
        const U u_en = _u_en(x);

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

    auto P1_dot(const _state_type &x)
    {
        const auto D = U(0.0); // enteral feed rate TODO: get from network
        return -d1 * x[fn_P1] + D;
    }

    auto P2_dot(const _state_type &x)
    {
        return -P_min(x) + d1 * x[fn_P1];
    }

    void copy(const ChaseIcing &other)
    {
        data = other.data;
        insulin_rate = other.insulin_rate;
        dextrose_rate = other.dextrose_rate;
    }

    /* data is the container which contains the most up to date representation of the model's state
     *
     * Should only be accessed using the enum Function data type
     */
    _state_type data;
public:
    using state_type = _state_type;
    void operator() (const state_type &x, state_type &dxdt, const U t)
    {
        dxdt[fn_G] = G_dot(x);
        dxdt[fn_Q] = Q_dot(x);
        dxdt[fn_I] = I_dot(x);
        dxdt[fn_P1] = P1_dot(x);
        dxdt[fn_P2] = P2_dot(x);
    }

    auto glucose()
    {
        return data[fn_G];
    }

    auto q()
    {
        return data[fn_Q];
    }

    auto i()
    {
        return data[fn_I];
    }

    auto p1()
    {
        return data[fn_P1];
    }

    auto p2()
    {
        return data[fn_P2];
    }

    int run(U time_start
            , U time_end
            , U dt
            , U insulin_rate_mUpermin
            , U dextrose_rate_mmolpermin)
    {

        insulin_rate = insulin_rate_mUpermin;
        dextrose_rate = dextrose_rate_mmolpermin;

        /* model doesnt modify inplace, so we make copies of its data
         * for integration */
        auto x = data;
        auto step = stepper;

        integrate_const(step, *this, x, time_start, time_end, dt);
        
        data = x;
        stepper = step;
        
        return 0;
    }

    /* run the model until time_end using default rates and time step
     */
    int run(U time_end)
    {
        return run(0.0, time_end, 0.1, insulin_rate, dextrose_rate);
    }

    /* run the model with implicitly provided rates
     */
    int run(U time_start, U time_end, U dt)
    {
        return run(time_start, time_end, dt, insulin_rate, dextrose_rate);
    }

    /* insulin_rate is the exogenous IV insulin rate, in mU/min (not mmol/min as listed in [1])
     *
     * Available for public modification, but use run() for most purposes.
     */
    U insulin_rate;

    /* dextrose_rate is the exogenous IV dextrose rate, in mmol/min
     *
     * Available for public modification, but use run() for most purposes.
     */
    U dextrose_rate;

    /* stepper is the Boost integration type, we use RK4 but you may use any you like
     */
    runge_kutta4<state_type> stepper;

    /* public constructor
     */
    ChaseIcing(state_type _data) : data(_data) {};

    /* copy constructor
     */
    ChaseIcing(const ChaseIcing &other)
    {
        copy(other);
    }

    /* assignment operator
     */
    ChaseIcing& operator=(const ChaseIcing &other)
    {
        copy(other);
        return *this;
    }
};

