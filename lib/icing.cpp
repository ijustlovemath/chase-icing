#include "ChaseIcing.hpp"


const double initial_glucose = 6.1;
using ModelType = ChaseIcing<double>;
static ModelType model;// = ModelType(initial_state);

extern int init_model(double initial_glucose)
{
	static ModelType::state_type initial_state = {initial_glucose, 19.0, 67.0, 277.54, 11.1};
	initial_state[0] = initial_glucose;
	model = (ModelType(initial_state));
	if(initial_glucose < 0.0 || initial_glucose > 25.0) {
	    return 3;
	}
	return 0;
}

extern int run_model(double time_start, double time_end, double dt, double insulin_rate, double dextrose_rate)
{
    return model.run(time_start, time_end, dt, insulin_rate, dextrose_rate);
}

extern int set_enteral(double *enteral)
{
	model.enteral_feed_rate = *enteral;
	return 0;
}

extern int latest_glucose(double *dest)
{
	*dest = model.glucose();
	return 0;
}
