#include <vector> /* std::vector */
#include <utility> /* std::pair */
#include <algorithm> /* std::lower_bound, std::sort */

#include "XYParser.hpp"

#include <boost/range/combine.hpp> /* boost::combine */
template<typename U = double>
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
