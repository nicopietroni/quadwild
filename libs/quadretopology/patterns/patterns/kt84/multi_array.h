#pragma once
#include <vector>
#include <Eigen/Core>

namespace kt84 {
    namespace internal {
        template <typename T, int N>
        struct multi_array_helper;
    }

template <typename T, int N>
struct multi_array {
    friend struct internal::multi_array_helper<T, N>;
    
    typedef Eigen::Matrix<int, N, 1> VectorNi;

    static int total_size(const VectorNi& size) {
        int result = 1;
        for (int i = 0; i < N; ++i) result *= size[i];
        return result;
    }
    
    // ctor()
    multi_array()
        : size_(VectorNi::Zero())
    {}
    // ctor(size)
    multi_array(const VectorNi& size__)
        : size_(size__)
        , data(total_size(size_))
    {}
    // ctor(size, val)
    multi_array(const VectorNi& size__, const T& val)
        : size_(size__)
        , data(total_size(size_), val)
    {}
    // ctor(sx, ...) specialized for 2D & 3D
    multi_array(int sx, int sy)         { internal::multi_array_helper<T, N>::ctor(sx, sy,     size_, data); }
    multi_array(int sx, int sy, int sz) { internal::multi_array_helper<T, N>::ctor(sx, sy, sz, size_, data); }
    multi_array(int sx, int sy,         const T& val) { internal::multi_array_helper<T, N>::ctor(sx, sy,     val, size_, data); }
    multi_array(int sx, int sy, int sz, const T& val) { internal::multi_array_helper<T, N>::ctor(sx, sy, sz, val, size_, data); }
    
    void clear() {
        size_.setZero();
        data.clear();
    }
    
    VectorNi size() const { return size_; }
    int size(int coord) const { return size_[coord]; }
    int total_size() const { return total_size(size_); }
    bool empty() const { return data.empty(); }
    
    // width() and height() for 2D & 3D, depth() for 3D
    int width () const { return internal::multi_array_helper<T, N>::width (*this); }
    int height() const { return internal::multi_array_helper<T, N>::height(*this); }
    int depth () const { return internal::multi_array_helper<T, N>::depth (*this); }
    
    int unroll_index(const VectorNi& index) const {
        int result = index[N - 1];
        for (int i = N - 2; i >= 0; --i) {
            result *= size_[i];
            result += index[i];
        }
        return result;
    }
    int unroll_index(int ix, int iy)         const { return internal::multi_array_helper<T, N>::unroll_index(*this, ix, iy); }
    int unroll_index(int ix, int iy, int iz) const { return internal::multi_array_helper<T, N>::unroll_index(*this, ix, iy, iz); }
    
    // operator(index)
    const T& operator()(const VectorNi& index) const { return data[unroll_index(index)]; }
          T& operator()(const VectorNi& index)       { return data[unroll_index(index)]; }
    // operator(ix, ...) specialized for 2D & 3D
    const T& operator()(int ix, int iy)         const { return internal::multi_array_helper<T, N>::at(*this, ix, iy); }
          T& operator()(int ix, int iy)               { return internal::multi_array_helper<T, N>::at(*this, ix, iy); }
    const T& operator()(int ix, int iy, int iz) const { return internal::multi_array_helper<T, N>::at(*this, ix, iy, iz); }
          T& operator()(int ix, int iy, int iz)       { return internal::multi_array_helper<T, N>::at(*this, ix, iy, iz); }
    // operator[int]
    const T& operator[](int i) const { return data[i]; }
          T& operator[](int i)       { return data[i]; }
    
    // resize(size)
    void resize(const VectorNi& size__) {
        size_ = size__;
        data.resize(total_size(size_));
    }
    // resize(size, val)
    void resize(const VectorNi& size__, const T& val) {
        size_ = size__;
        data.resize(total_size(size_), val);
    }
    // resize(sx, ...) specialized for 2D & 3D
    void resize(int sx, int sy)         { internal::multi_array_helper<T, N>::resize(*this, sx, sy); }
    void resize(int sx, int sy, int sz) { internal::multi_array_helper<T, N>::resize(*this, sx, sy, sz); }
    void resize(int sx, int sy,         const T& val) { internal::multi_array_helper<T, N>::resize(*this, sx, sy, val); }
    void resize(int sx, int sy, int sz, const T& val) { internal::multi_array_helper<T, N>::resize(*this, sx, sy, sz, val); }
    
    // begin/end
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typename iterator       begin()       { return data.begin(); }
    typename iterator       end  ()       { return data.end  (); }
    typename const_iterator begin() const { return data.begin(); }
    typename const_iterator end  () const { return data.end  (); }
    typename const_iterator cbegin() const { return begin(); }
    typename const_iterator cend  () const { return end  (); }
    
private:
    VectorNi size_;
    std::vector<T> data;
};

    namespace internal {
        template <typename T>
        struct multi_array_helper<T, 2> {
            static void ctor(int sx, int sy, Eigen::Vector2i& size, std::vector<T>& data) {
                size = Eigen::Vector2i(sx, sy);
                data.resize(multi_array<T, 2>::total_size(size));
            }
            static void ctor(int sx, int sy, const T& val, Eigen::Vector2i& size, std::vector<T>& data) {
                size = Eigen::Vector2i(sx, sy);
                data.resize(multi_array<T, 2>::total_size(size), val);
            }
            static int width (const multi_array<T, 2>& m) { return m.size_.x(); }
            static int height(const multi_array<T, 2>& m) { return m.size_.y(); }
            static int unroll_index(const multi_array<T, 2>& m, int ix, int iy) { return m.unroll_index(Eigen::Vector2i(ix, iy)); }
            static const T& at(const multi_array<T, 2>& m, int ix, int iy) { return m(Eigen::Vector2i(ix, iy)); }
            static       T& at(      multi_array<T, 2>& m, int ix, int iy) { return m(Eigen::Vector2i(ix, iy)); }
            static void resize(multi_array<T, 2>& m, int sx, int sy)               { m.resize(Eigen::Vector2i(sx, sy)); }
            static void resize(multi_array<T, 2>& m, int sx, int sy, const T& val) { m.resize(Eigen::Vector2i(sx, sy), val); }
        };
        template <typename T>
        struct multi_array_helper<T, 3> {
            static void ctor(int sx, int sy, int sz, Eigen::Vector3i& size, std::vector<T>& data) {
                size = Eigen::Vector3i(sx, sy, sz);
                data.resize(multi_array<T, 3>::total_size(size));
            }
            static void ctor(int sx, int sy, int sz, const T& val, Eigen::Vector3i& size, std::vector<T>& data) {
                size = Eigen::Vector3i(sx, sy, sz);
                data.resize(multi_array<T, 3>::total_size(size), val);
            }
            static int width (const multi_array<T, 3>& m) { return m.size_.x(); }
            static int height(const multi_array<T, 3>& m) { return m.size_.y(); }
            static int depth (const multi_array<T, 3>& m) { return m.size_.z(); }
            static int unroll_index(const multi_array<T, 3>& m, int ix, int iy, int iz) { return m.unroll_index(Eigen::Vector3i(ix, iy, iz)); }
            static const T& at(const multi_array<T, 3>& m, int ix, int iy, int iz) { return m(Eigen::Vector3i(ix, iy, iz)); }
            static       T& at(      multi_array<T, 3>& m, int ix, int iy, int iz) { return m(Eigen::Vector3i(ix, iy, iz)); }
            static void resize(multi_array<T, 3>& m, int sx, int sy, int sz)               { m.resize(Eigen::Vector3i(sx, sy, sz)); }
            static void resize(multi_array<T, 3>& m, int sx, int sy, int sz, const T& val) { m.resize(Eigen::Vector2i(sx, sy, sz), val); }
        };
    }
}
