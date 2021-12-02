#pragma once

namespace kt84 {

#ifdef NDEBUG
    inline void myassert(bool test) {}
#else
    inline void myassert(bool test) { test ? 1 : (*reinterpret_cast<int*>(0) = 1); }
#endif

}
