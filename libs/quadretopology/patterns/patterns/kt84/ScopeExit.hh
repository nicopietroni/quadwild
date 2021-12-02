#pragma once

namespace kt84 {
    template <typename Func>
    struct ScopeExit {
        ScopeExit(Func func) : func(func) {}
        ~ScopeExit() { func(); }
        Func func;
    };
    template <typename Func>
    inline ScopeExit<Func> make_ScopeExit(Func func) { return ScopeExit<Func>(func); }
}
