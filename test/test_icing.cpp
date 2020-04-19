#include <vector> /* representing state */
#include <algorithm> /* for numeric min/max */
#include <cmath> /* log, exp, fmod */
#include <iostream> /* cout */
#include <string> /* string */
#include <ios> /* std::setw */

#include <imt/imt.h>

#include "ChaseIcing.hpp"
#include "InverseCDFProcess.hpp"
#include "Interpolator.hpp"

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

int run_model(Chase &model
        , double time_start
        , imt_float_t time_step
        , imt_float_t ins
        , imt_float_t dex)
{

    const double dt = 0.1;
    return model.run(time_start, time_start + time_step, dt, ins, dex);
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
    Chase::state_type x = {initial_glucose, 19.0, 67.0, 277.54, 11.1};
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
        if(err) {
            std::cout << imt_strerror(err);
            running = 0;
        }

        if(enable_controller) {
            /* get insulin rate */
            err = imt_insulin_prescription_in(ctx, INS_units
                    , rates_normalized_by_weight
                    , &insulin_rate
            );
            if(err) {
                std::cout << imt_strerror(err);
                running = 0;
            }

            /* get dextrose rate */
            err = imt_dextrose_prescription_in(ctx, DEX_units
                    , rates_normalized_by_weight
                    , &dextrose_rate
            );
            if(err) {
                std::cout << imt_strerror(err);
                running = 0;
            }
        }

        /* get time step, in minutes */
        err = imt_cycle_length_min(ctx, &time_step);
        if(err) {
            std::cout << imt_strerror(err);
            running = 0;
        }

        /* log the current data to stdout */
        double ml_ins = -1.0, ml_dex = -1.0;
        imt_insulin_prescription_in(ctx, "mL/hr", 0, &ml_ins);
        imt_dextrose_prescription_in(ctx, "mL/hr", 0, &ml_dex);
        print_row(time, model, ml_ins, ml_dex);

        /* run the model you've chosen for one cycle */
        err = run_model(model, time, time_step, insulin_rate, dextrose_rate);
        if(err) {
            std::cout << imt_strerror(err);
            running = 0;
        }

        time += time_step;

        /* get the latest glucose value from the model */
        /* stubbed for now */
        imt_float_t next_glucose = imt_mmolperL_to_mgperdL(model.glucose());

        /* run reassignment before running the controller again */
        err = imt_exec_reassignment(ctx, next_glucose);
        if(err) {
            std::cout << imt_strerror(err);
            running = 0;
        }

        if(running)
            running = time <= tend;

    }
}

void test_interpolator_inversecdf(void)
{
    auto interp = Interpolator<double>("../data/ICING BE SI Data.csv");
    interp.set_bounds(std::make_pair<double, double>(0.0, 1.0));

    auto generator = InverseCDFProcess<double>("../data/ICING BE SI Data.csv");

    for(int i = 0; i < 30; ++i) {
        std::cout << generator.generate() << std::endl;
    }
}

int main(void)
{
    run_controller_with_model(true);
    return 0;
}
