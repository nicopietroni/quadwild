#pragma once
#include <cmath>
#include <cstdarg>
#include <limits>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <iostream>
#include <random>
#include <ctime>

namespace kt84 {

namespace util {
    template <typename T> inline T lexical_cast(const std::string& str) {
        T ret;
        std::stringstream ss;
        ss.str(str);
        ss >> ret;
        return ret;
    }
    template <typename T> inline void append_str(std::string& str, const T& value) {
        std::stringstream ss;
        ss << str << value;
        str = ss.str();
    }
    template <typename T> inline T& incr_mod(T& value, const T& mod) { return value = (value + 1) % mod; }
    template <typename T> inline T& decr_mod(T& value, const T& mod) { return value = (value + mod - 1) % mod; }
    inline bool ensure_filename_extension(std::string& filename, const std::string& extension_including_period) {
        auto n = extension_including_period.size();
        if (!boost::iequals(filename.substr(filename.size() - n, n), extension_including_period)) {
            filename += extension_including_period;
            return false;
        }
        return true;
    }
    inline void flip_bool(bool& b) { b = !b; }
    template <typename T> inline void update_if_bigger (T& val_current, const T& val_new) { if (val_current < val_new) val_current = val_new; }
    template <typename T> inline void update_if_smaller(T& val_current, const T& val_new) { if (val_current > val_new) val_current = val_new; }
    template <typename T> inline T clamp(const T& value, const T& range_min, const T& range_max) { return std::max<T>(std::min<T>(value, range_max), range_min); }
    inline double asin_clamped(double sine)   { return std::asin(clamp(sine  , -1.0, 1.0)); };
    inline double acos_clamped(double cosine) { return std::acos(clamp(cosine, -1.0, 1.0)); };
    inline double pi() { return std::acos(-1); }
    inline double dbl_max() { return std::numeric_limits<double>::max(); }
    inline int    int_max() { return std::numeric_limits<int   >::max(); }
    template <typename T> inline T squared(const T& value) { return value * value; }
    template <typename T> inline T cubed(const T& value) { return value * value * value; }
    template <typename T> inline T normalize(const T& value, const T& value_min, const T& value_max) { return (value - value_min) / (value_max - value_min); }
    
    // max/min of more than two inputs. this is already possible in C++11 using initializer lists, but unfortunately it's still unsupported in VC11.
    template <typename T> inline T max(const T& arg0, const T& arg1, const T& arg2) { return std::max<T>(std::max<T>(arg0, arg1), arg2); }
    template <typename T> inline T min(const T& arg0, const T& arg1, const T& arg2) { return std::min<T>(std::min<T>(arg0, arg1), arg2); }
    template <typename T> inline T max(const T& arg0, const T& arg1, const T& arg2, const T& arg3) { return std::max<T>(max(arg0, arg1, arg2), arg3); }
    template <typename T> inline T min(const T& arg0, const T& arg1, const T& arg2, const T& arg3) { return std::min<T>(min(arg0, arg1, arg2), arg3); }
    template <typename T> inline T max(const T& arg0, const T& arg1, const T& arg2, const T& arg3, const T& arg4) { return std::max<T>(max(arg0, arg1, arg2, arg3), arg4); }
    template <typename T> inline T min(const T& arg0, const T& arg1, const T& arg2, const T& arg3, const T& arg4) { return std::min<T>(min(arg0, arg1, arg2, arg3), arg4); }

    // not so great with va_arg, but just for convenience...
    template <typename T>
    inline std::vector<T> make_vector(int size, ...) {
        va_list args;
        va_start(args, size);
        std::vector<T> result;
        result.reserve(size);
        for (int i = 0; i < size; ++i)
            result.push_back(va_arg(args, T));
        va_end(args);
        return result;
    }
    template <typename T>
    inline std::list<T> make_list(int size, ...) {
        va_list args;
        va_start(args, size);
        std::list<T> result;
        for (int i = 0; i < size; ++i)
            result.push_back(va_arg(args, T));
        va_end(args);
        return result;
    }
    template <typename T>
    inline T simple_prompt(const std::string& name) {
        T value;
        std::cout << name << ":";
        std::cin >> value;
        return value;
    }
    template <typename T>
    inline T simple_prompt_array(const std::string& name, int n) {
        T value;
        std::cout << name << ":";
        for (int i = 0; i < n; ++i)
            std::cin >> value[i];
        return value;
    }
    template <typename T>
    inline T simple_prompt(const std::string& name, const T& current_value) {
        T value;
        std::cout << name << " (current=" << current_value << "):";
        std::cin >> value;
        return value;
    }
    template <typename T>
    inline void swap_pair(std::pair<T, T>& p) { std::swap(p.first, p.second); }
    template <typename T>
    inline bool is_pair_equal(const std::pair<T, T>& p1, const std::pair<T, T>& p2) {       // compare p1 and p2 while ignoring ordering
        return p1.first == p2.first  && p1.second == p2.second ||
               p1.first == p2.second && p1.second == p2.first;
    }
    inline std::mt19937& random_engine() {
        static std::mt19937 engine(static_cast<unsigned int>(std::time(nullptr)));
        return engine;
    }
    template <typename T>
    T random_real(T range_min, T range_max) {   // range: [min, max)
        static std::uniform_real_distribution<T> dist(0, 1);       // range: [0, 1)
        T t = dist(random_engine());
        return (1 - t) * range_min + t * range_max;
    }
    inline double random_double(double range_min, double range_max) { return random_real<double>(range_min, range_max); }
    inline float  random_float (float  range_min, float  range_max) { return random_real<float >(range_min, range_max); }
    template <typename T>
    T random_integral(T range_min, T range_max) {               // range: [min, max]
        return range_min + static_cast<T>(random_double(0, 1) * (range_max - range_min + 1));
    }
    inline int  random_int (int  range_min, int  range_max ) { return random_integral<int > (range_min, range_max); }
    inline int  random_long(long range_min, long range_max ) { return random_integral<long> (range_min, range_max); }
    inline void srand(unsigned long seed) {
        random_engine().seed(seed);
    }
    template <typename T1, typename T2>
    inline void unfold_pair(const std::pair<T1, T2>& p, T1& val1, T2& val2) { val1 = p.first; val2 = p.second; }
    template <typename TEnum, typename TInt>
    inline TEnum add_enum(const TEnum& e, TInt n) { return static_cast<TEnum>(static_cast<TInt>(e) + n); }
    template <typename Value, typename IStream>
    Value query_istream(IStream& stream) {
        Value value;
        stream >> value;
        return value;
    }
}

}
