#include <random>
#include <vector>
#include <string>
#include "Interpolator.hpp"
template<typename U = double>
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
