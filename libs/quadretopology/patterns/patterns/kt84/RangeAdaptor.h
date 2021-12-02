#pragma once

namespace kt84 {
    template <typename TIterator>
    struct RangeAdaptor {
        RangeAdaptor(TIterator begin_, TIterator end_) : iter_begin(begin_), iter_end(end_) {}
        TIterator begin() const { return iter_begin; }
        TIterator end  () const { return iter_end  ; }
    private:
        TIterator iter_begin, iter_end;
    };
    template <typename TIterator>
    inline RangeAdaptor<TIterator> make_RangeAdaptor(TIterator begin_, TIterator end_) { return RangeAdaptor<TIterator>(begin_, end_); }
}
