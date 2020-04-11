#include <vector> /* representing state */
#include <algorithm> /* for numeric min/max */
#include <cmath> /* exp, fmod */
#include <iostream> /* cout */
#include <string> /* string */

#include <boost/numeric/odeint.hpp>
#include <ios>
using namespace boost::numeric::odeint;

#include <imt/imt.h>
#include <imt/libimt.h>

typedef double icing_float;
//typedef std::vector<icing_float> state_type;

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

    /* a bunch of constants that should be in the namespace */
    const U n_I = 0.0075; // 1/min, Table II, [1]
    const U d1 = -std::log(0.5)/20.0; // 1/min, Table II, [1]
    const U d2 = -std::log(0.5)/100.0; // 1/min, Table II, [1]
    const U P_max = 6.11; // mmol/min, Table II, [1]

    auto P_min(const std::vector<U> &x)
    {
        return std::min(d2 * x[fn_P2], P_max);
    }

    auto alpha_decay(const U variable, const U decay_parameter)
    {
        return variable / (1.0 + decay_parameter * variable);
    }

    auto Q_frac(const std::vector<U> &x)
    {
        const U a_G = 1.0 / 65.0; // Table II, [1]
        const U Q = x[fn_Q];
        return alpha_decay(Q, a_G);
    }

    // PN(t) -> Parenteral nutrition input, eg IV dextrose
    auto _P(const std::vector<U> &x)
    {
        const U _PN_ext = dextrose_rate; // TODO -> derive this over the network

        return P_min(x) + _PN_ext;
    }

    auto _u_en(const std::vector<U> &x)
    {
        const U k1 = 14.9; // mU * l / mmol/min, Table II, [1]
        const U k2 = -49.9; // mU/min, Table II, [1]
        const U u_min = 16.7; // mU/min, Table II, [1]
        const U u_max = 266.7; // mU/min, Table II, [1]
        const U G = x[fn_G];
        return std::min(std::max(u_min, k1 * G + k2), u_max);
    }

    auto G_dot(const std::vector<U> &x)
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

    auto Q_dot(const std::vector<U> &x)
    {
        const U I = x[fn_I];
        const U Q = x[fn_Q];
        const U n_C = 0.0075; // 1/min, Table II, [1]

        return n_I * (I - Q) - n_C * Q_frac(x);
    }

    auto I_dot(const std::vector<U> &x)
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

    auto P1_dot(const std::vector<U> &x)
    {
        const auto D = U(0.0); // enteral feed rate TODO: get from network
        return -d1 * x[fn_P1] + D;
    }

    auto P2_dot(const std::vector<U> &x)
    {
        return -P_min(x) + d1 * x[fn_P1];
    }

public:
    using state_type = std::vector<U>;
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

    double insulin_rate;
    double dextrose_rate;
    runge_kutta4<state_type> stepper;
    state_type data;

    ChaseIcing(state_type _data) : data(_data) {};
};

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
/* largely taken from https://stackoverflow.com/a/415553/2289030 */
template<typename U = icing_float>
void load_xy_data(const std::string &filename, std::vector<U> &x, std::vector<U> &y)
{
    std::ifstream data(filename);

    std::string line;
    while(std::getline(data, line))
    {
        std::stringstream linestream(line);
        std::string cell;
        int i = 0;
        while(std::getline(linestream, cell, ','))
        {
            const U element = static_cast<U>(std::stod(cell));
            switch(i) {
            case 0:
                x.push_back(element);
                break;
            case 1:
                y.push_back(element);
                break;
            default: ;
            }
            ++i;
        }
        if(i < 2) { 
            // avoid mismatch if line is malformated
            x.pop_back();
            std::cerr << "problem loading data from " << filename << std::endl;
        }
    }
}

#include <vector> /* std::vector */
#include <utility> /* std::pair */
#include <algorithm> /* std::lower_bound, std::sort */
#include <boost/range/combine.hpp> /* boost::combine */
template<typename U = icing_float>
class Interpolator {
    std::vector<std::pair<U, U>> _table;
    const U INF = static_cast<U>(1e100);

    U _lower_bound = -INF;
    U _upper_bound = INF;

    void init(const std::vector<U> &x, const std::vector<U> &y)
    {
        // taken from here: https://stackoverflow.com/a/8513803/2289030
        for(auto tup : boost::combine(x, y)) {
            U x, y;
            boost::tie(x, y) = tup;

            _table.push_back(std::make_pair(x, y));
        }

        std::sort(_table.begin(), _table.end());
    }

    void init(const std::vector<U> &x, const std::vector<U> &y, std::pair<U, U> bounds)
    {
        init(x, y);
        set_bounds(bounds);
    }
    
    void copy(const Interpolator &other)
    {
        _table = other._table;
        _lower_bound = other._lower_bound;
        _upper_bound = other._upper_bound;    
    }
public:
    U get(U x)
    {
        // taken from here: https://stackoverflow.com/a/11675205/2289030
        // Assumes that "table" is sorted by .first

        // Check if x is out of bound
        if (x > _table.back().first) return _upper_bound;
        if (x < _table[0].first) return _lower_bound;

        typename std::vector<std::pair<U, U>>::iterator it, it2;
        // INFINITY is defined in math.h in the glibc implementation
        it = std::lower_bound(_table.begin(), _table.end(), std::make_pair(x, -INF));
        // Corner case
        if (it == _table.begin()) return it->second;
        it2 = it;
        --it2;
        return it2->second + ((it->second - it2->second)
                                *(x - it2->first)/(it->first - it2->first));
    }

    void set_bounds(const std::pair<U, U> &bounds)
    {
        _lower_bound = bounds.first;
        _upper_bound = bounds.second;
    }

    U interpolate(U x)
    {
        return get(x);
    }

    Interpolator(const std::vector<U> &x, const std::vector<U> &y)
    {
        init(x, y);
    }

    Interpolator(const std::vector<U> &x, const std::vector<U> &y
            , std::pair<U,U> bounds)
    {
        init(x, y, bounds);
    }

    Interpolator(const std::string &filename)
    {
        std::vector<U> x, y;
        load_xy_data<U>(filename, x, y);
        
        init(x, y);
    }

    Interpolator() {};

    Interpolator(const Interpolator &other)
    {
        copy(other);
    }

    Interpolator& operator=(const Interpolator &rhs)
    {
        copy(rhs);
        return *this;
    }

};

#include <random>
//#include "Interpolator.hpp"
template<typename U = icing_float>
class InverseCDFProcess {

    std::vector<U> _x;
    std::vector<U> _y;

    std::random_device _r;
    std::seed_seq seed{_r(), _r(), _r(), _r(), _r(), _r(), _r()};
    //std::default_random_engine _generator;
    std::mt19937 _generator{seed};
    std::uniform_real_distribution<U> _distribution{0.0, 1.0};

    Interpolator<U> _interp;

    const std::pair<U,U> bounds = std::make_pair<U,U>(0.0, 1.0);
    
    void copy(const InverseCDFProcess &other)
    {
        _interp = other._interp;
        _x = other._x;
        _y = other._y;
    }

public:
    InverseCDFProcess(const std::string &filename)
    {
        load_xy_data<U>(filename, _x, _y);
        _interp = Interpolator<U>(_y, _x, bounds);
    }
    InverseCDFProcess(const std::vector<U> &x, const std::vector<U> &y) : 
        _x(x), _y(y)
        , _interp(Interpolator<U>(y, x, bounds)) {};

    InverseCDFProcess(const InverseCDFProcess &other)
    {
        copy(other);
    }

    InverseCDFProcess& operator=(const InverseCDFProcess &other)
    {
        copy(other);
        return *this;
    }

    U generate()
    {
        // https://www.cplusplus.com/reference/random/uniform_real_distribution/
        return _interp.interpolate(_distribution(_generator));
    }
};

using Chase = ChaseIcing<>;

void print_row(double t, Chase model, double i, double d)
{
    static bool first_call = true;
    const auto width = 25;
    const auto g = model.glucose();
    const auto q = model.q();
    const auto I = model.i();
    const auto p1 = model.p1();
    const auto p2 = model.p2();


    if(first_call) {
        std::cout << std::setw(width) << "Time (min)";
        std::cout << std::setw(width) << "Glucose (mg/dL)";
        std::cout << std::setw(width) << "Insulin (mL/hr)";
        std::cout << std::setw(width) << "Dextrose (mL/hr)";
        std::cout << std::setw(width) << "Q";
        std::cout << std::setw(width) << "I";
        std::cout << std::setw(width) << "P1";
        std::cout << std::setw(width) << "P2";
        std::cout << std::endl;
        first_call = false;
    }
    std::cout << std::setw(width) << t;
    std::cout << std::setw(width) << imt_mmolperL_to_mgperdL(g);
    std::cout << std::setw(width) << i;
    std::cout << std::setw(width) << d;
    std::cout << std::setw(width) << q;
    std::cout << std::setw(width) << I;
    std::cout << std::setw(width) << p1;
    std::cout << std::setw(width) << p2;
    std::cout << std::endl;

}

using vec5 = Chase::state_type;

int run_model(Chase &model
        , double time_start
        , imt_float_t time_step
        , imt_float_t ins
        , imt_float_t dex)
{

    const double dt = 0.1;
    model.insulin_rate = ins;
    model.dextrose_rate = dex;
    /* model doesnt modify inplace, so we make copies of its data
     * for integration */
    auto x = model.data;
    auto step = model.stepper;

    integrate_const(step, model, x, time_start, time_start+time_step, dt);
    
    model.data = x;
    model.stepper = step;
    
    return 0;
}

void run_controller_with_model(bool enable_controller)
{
    const double tend = 1e5;
	
    const int rates_normalized_by_weight = 0;

    const imt_float_t weight = 80.0; /* kg */
    const imt_float_t INS_conc = 1.0; /* 1 U/mL */
    const imt_float_t DEX_conc = 20.0; /* D20% */
    const imt_float_t control_min = 4.4; /* mmol/L */
    const imt_float_t control_max = 9.0; /* mmol/L */
    const imt_float_t initial_glucose = 6.1; /* mmol/L */

    const char * INS_units = "mU/min";
    const char * DEX_units = "mmol/min";

    imt_context_t context, *ctx = &context;
    double insulin_rate=0.0, dextrose_rate=0.0;
    double time_step;

    /* initialize model */
    /* plasma insulin 19 +/- 2 mU/L */
    /* interstitial insulin 67 +/- 19 mU/L */
    /* glucose in stomach 50g postprandial = 277.54 mmol */
    /* glucose in gut ? maybe 2g = 11.1 mmol */
    vec5 x = {initial_glucose, 19.0, 67.0, 277.54, 11.1};
    Chase model(x);  

    /* sets up controller to control to a range of 4.4-9.0 mmol/L */
    imt_exec_init(ctx 
        , imt_mmolperL_to_mgperdL(initial_glucose)
        , weight
        , imt_mmolperL_to_mgperdL(control_min)
        , imt_mmolperL_to_mgperdL(control_max)
        , INS_conc
        , DEX_conc
    );

    int running = 1, err = I_OKAY;
    double time = 0.0;

    /* proper error handling is NOT done here */
    while(running) {
        /* run the controller for one cycle */
        err = imt_exec_controller(ctx);
        if(err)
            running = 0;

        if(enable_controller) {
            /* get insulin rate */
            err = imt_insulin_prescription_in(ctx, INS_units
                    , rates_normalized_by_weight
                    , &insulin_rate
            );
            if(err)
                running = 0;

            /* get dextrose rate */
            err = imt_dextrose_prescription_in(ctx, DEX_units
                    , rates_normalized_by_weight
                    , &dextrose_rate
            );
            if(err)
                running = 0;
        }

        /* get time step, in minutes */
        err = imt_cycle_length_min(ctx, &time_step);
        if(err)
            running = 0;

        /* log the current data to stdout */
        double ml_ins = -1.0, ml_dex = -1.0;
        imt_insulin_prescription_in(ctx, "mL/hr", 0, &ml_ins);
        imt_dextrose_prescription_in(ctx, "mL/hr", 0, &ml_dex);
        print_row(time, model, ml_ins, ml_dex);

        /* run the model you've chosen for one cycle */
        err = run_model(model, time, time_step, insulin_rate, dextrose_rate);
        if(err)
            running = 0;

        time += time_step;

        /* get the latest glucose value from the model */
        /* stubbed for now */
        imt_float_t next_glucose = imt_mmolperL_to_mgperdL(model.glucose());

        /* run reassignment before running the controller again */
        err = imt_exec_reassignment(ctx, next_glucose);
        if(err)
            running = 0;

        if(running)
            running = time <= tend;

    }
}

int main(void)
{

    auto interp = Interpolator<double>("../data/ICING BE SI Data.csv");
    interp.set_bounds(std::make_pair<double, double>(0.0, 1.0));

    std::cout << interp.interpolate(0.5);
    auto generator = InverseCDFProcess<double>("../data/ICING BE SI Data.csv");

    for(int i = 0; i < 30; ++i) {
        std::cout << generator.generate() << std::endl;
    }
    return 0;
}
