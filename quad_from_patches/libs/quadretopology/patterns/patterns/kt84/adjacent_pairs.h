#pragma once
#include "RangeAdaptor.h"

namespace kt84 {
    namespace internal {
        template<typename Container, typename Iterator, typename Element>
        struct AdjacentPairIter {
            AdjacentPairIter(Container& parent_, Iterator pos_)
                : parent(&parent_), pos(pos_)
            {}
            // minimal required member functions
            typedef std::pair<Element&, Element&> DerefType;
            DerefType operator*() const {
                Iterator pos_next(pos);
                auto& first  = *pos;
                auto& second = *(++pos_next == parent->end() ? parent->begin() : pos_next);
                return DerefType(first, second);
            }
            AdjacentPairIter& operator++() {
                ++pos;
                return *this;
            }
            bool operator!=(const AdjacentPairIter& rhs) const { return parent != rhs.parent || pos != rhs.pos; }
        private:
            Container* parent;
            Iterator pos;
        };
    }
    template<class Container>
    inline auto adjacent_pairs(Container& c, bool is_loop) -> RangeAdaptor<internal::AdjacentPairIter<Container, decltype(c.begin()), decltype(*c.begin())>> {
        typedef decltype( c.begin()) Iterator;
        typedef decltype(*c.begin()) Element;
        internal::AdjacentPairIter<Container, Iterator, Element> iter_begin(c, c.begin());
        internal::AdjacentPairIter<Container, Iterator, Element> iter_end  (c, is_loop ? c.end() : --c.end());
        return make_RangeAdaptor(iter_begin, iter_end);
    }
}
