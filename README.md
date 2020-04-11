# chase-icing
C++ Implementation of the 2016 Chase ICING model

# Building

This is a CMake project. To build the test executable, unfortunately you'll need access to a private repository.

If you did have access to the repository, the build would be as simple as:

`mkdir build && cd build && cmake .. && make` (on linux)

# Usage

The most important aspect of this repo is not the test executable (which shows how to attach a controller to the model), but the model itself.

The model is a header-only library, contained in `include/ChaseIcing.hpp`. This contains a templated ChaseIcing<U> class.

It can be instantiated in a default format using:

```
ChaseIcing<>::state_type initial_state = {...};
auto model = ChaseIcing<>(initial_state);
```

The model can be run for a specified number of minutes `([time_start, time_end])`, with a specified update interval `dt`, and exogenous IV insulin flow `insulin_rate` (in mU/min) and exogenous IV dextrose flow  `dextrose_rate` (in mmol/min) using `ChaseIcing<>().run()`:

`model.run(time_start, time_end, dt, insulin_rate, dextrose_rate)`

The resulting glucose, in mmol/L, can be obtained using `ChaseIcing<>().glucose()`:

`auto resulting_glucose = model.glucose()`

An example library, with a (mostly correct) initial state vector, can be found in `lib/icing.cpp`

# Future Work

In the future I plan to extend the model using randomly generated insulin sensitivity as derived from Fig 4 and Fig 5 of [1], to simulate patients properly from the original cohort of ICU patients. This will be implemented using the existing `InverseCDFProcess<U>` template class that's in the `include` folder. 

# References

[1] - TBME-01160-2016 
[2] - computer methods and programs in biomedicine 102 (2011) 192–205
