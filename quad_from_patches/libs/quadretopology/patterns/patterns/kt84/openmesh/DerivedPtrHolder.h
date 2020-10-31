#pragma once

namespace kt84 {
    template <typename TDerived, class TBase>
    struct DerivedPtrHolder {
        DerivedPtrHolder() : derived_ptr(static_cast<TDerived*>(this)) {}
        DerivedPtrHolder(const DerivedPtrHolder& rhs) : derived_ptr(static_cast<TDerived*>(this)) {}
        DerivedPtrHolder& operator=(const DerivedPtrHolder& rhs) {
            derived_ptr = static_cast<TDerived*>(this);
            return *this;
        }
        TDerived* derived_ptr;
    };
}
