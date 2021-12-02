#pragma once

namespace kt84 {

template <bool SumToOne>
class BaryCoordT;

template <>
class BaryCoordT<true> {
    double u;
    double v;
    
public:
    BaryCoordT()
        : u(1.0 / 3.0)
        , v(1.0 / 3.0)
    {}
    BaryCoordT(double w0, double w1, double w2) { set(w0, w1, w2); }
    
    void set(double w0, double w1, double w2) {
        double error = w0 + w1 + w2 - 1;
        assert (-error_tolerance() < error && error < error_tolerance());
        u = w1;
        v = w2;
    }
    inline double operator[](int index) const {
        return
            index == 1 ? u :
            index == 2 ? v :
            (1 - u - v);
    }
    bool is_all_positive() const { return u >= 0 && v >= 0 && (1 - u - v) >= 0; }
    
    template <typename T>
    T blend(const T& v0, const T& v1, const T& v2) const {
        return static_cast<T>(
            (1 - u - v) * v0 +
            u           * v1 +
            v           * v2);
    }
private:
    static double error_tolerance() { return 0.001; }
};

template <>
class BaryCoordT<false> {
    double u;
    double v;
    
public:
    BaryCoordT()
        : u(1.0 / 3.0)
        , v(1.0 / 3.0)
    {}
    BaryCoordT(double w0, double w1, double w2) { set(w0, w1, w2); }
    
    void set(double w0, double w1, double w2) {
        double error = w0 + w1 + w2;
        assert (-error_tolerance() < error && error < error_tolerance());
        u = w1;
        v = w2;
    }
    inline double operator[](int index) const {
        return
            index == 1 ? u :
            index == 2 ? v :
            (- u - v);
    }
    
    template <typename T>
    T blend(const T& v0, const T& v1, const T& v2) const {
        return static_cast<T>(
            (1 - u - v) * v0 +
            u           * v1 +
            v           * v2);
    }
private:
    static double error_tolerance() { return 0.001; }
};

typedef BaryCoordT<true > BaryCoord;
typedef BaryCoordT<false> BaryCoordZero;


}
