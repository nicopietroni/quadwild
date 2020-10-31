#pragma once

namespace kt84 {

template <typename T, typename DeleteFunc>
class AutoDeleter {
public:
    AutoDeleter(T* ptr, DeleteFunc delete_func) : ptr(ptr), delete_func(delete_func) {}
    AutoDeleter(const AutoDeleter& src) : ptr(src.ptr), delete_func(src.delete_func) {
        src.ptr = nullptr;
    }
    AutoDeleter& operator=(const AutoDeleter& src) {
        ptr         = src.ptr;
        delete_func = src.delete_func;
        src.ptr = nullptr;
        return *this;
    }
    ~AutoDeleter() {
        if (ptr != nullptr)
            delete_func(ptr);
    }
protected:
    T* ptr;
    DeleteFunc delete_func;
};
template <typename T, typename DeleteFunc>
inline AutoDeleter<T, DeleteFunc> make_AutoDeleter(T* ptr, DeleteFunc delete_func) {
    return AutoDeleter<T, DeleteFunc>(ptr, delete_func);
}

}
