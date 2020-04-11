#include "ChaseIcing.hpp"

const double initial_glucose = 6.1;
static const ChaseIcing<double>::state_type initial_state = {initial_glucose, 19.0, 67.0, 277.54, 11.1};
static auto model = ChaseIcing<double>(initial_state);

extern int run_model(double time_start, double time_end, double dt, double insulin_rate, double dextrose_rate)
{
    return model.run(time_start, time_end, dt, insulin_rate, dextrose_rate);
}
