#pragma once
#include <algorithm>
#include <numeric>
#include <boost/optional.hpp>
#include "RangeAdaptor.h"

namespace kt84 {
namespace container_util {
    template<class Map>
    inline boost::optional<typename Map::mapped_type> at_optional(Map& map, const typename Map::key_type& key) {
        auto p = map.find(key);
        if (p == map.end())
            return boost::none;
        else
            return p->second;
    }
    template<class Container>
    inline void unique(Container& c) {
        c.erase(std::unique(c.begin(), c.end()), c.end());
    }
    template<class Container>
    inline void reverse(Container& c) {
        std::reverse(c.begin(), c.end());
    }
    template<class Container>
    inline int mod_index(Container& c, int index) {
        while (index < 0) index += c.size();
        return index % c.size();
    }
    template<class Container>
    inline auto at_mod(Container& c, int index) -> decltype(c[index]) {
        return c[mod_index(c,index)];
    }
    template<class Container, class Predicate>
    inline void remove_if(Container& c, Predicate pred) {
        c.erase(std::remove_if(c.begin(), c.end(), pred), c.end());
    }
    template<class Container, class Type> 
    inline void bring_front(Container& c, const Type& val) {
        std::rotate(c.begin(), std::find(c.begin(), c.end(), val), c.end());
    }
    template<class Container, class Type> 
    inline auto find(Container& c, const Type& val) -> decltype(std::find(c.begin(), c.end(), val)) {
        return std::find(c.begin(), c.end(), val);
    }
    template<class Container, class Predicate> 
    inline auto find_if(Container& c, Predicate pred) -> decltype(std::find_if(c.begin(), c.end(), pred)) {
        return std::find_if(c.begin(), c.end(), pred);
    }
    template<class Container, class Type> 
    inline auto count(Container& c, const Type& val ) -> decltype(std::count(c.begin(), c.end(), val)) {
        return std::count(c.begin(), c.end(), val);
    }
    template<class Container, class Predicate> 
    inline auto count_if(Container& c, Predicate pred) -> decltype(std::count_if(c.begin(), c.end(), pred)) {
        return std::count_if(c.begin(), c.end(), pred);
    }
    template<class Container> 
    inline void sort(Container& c) {
        return std::sort(c.begin(), c.end());
    }
    template<class Container, class Predicate> 
    inline void sort(Container& c, Predicate pred) {
        return std::sort(c.begin(), c.end(), pred);
    }
    template<class Container> 
    inline auto erase_at(Container& c, int index) -> decltype(c.erase(c.begin())) {
        auto pos = c.begin();
        std::advance(pos, index);
        return c.erase(pos);
    }
    template<class Container, typename Element>
    inline Element accumulate(Container& c, const Element& val) {
        return std::accumulate(c.begin(), c.end(), val);
    }
    template<class Container, typename Element, typename BinaryOp>
    inline Element accumulate(Container& c, const Element& val, BinaryOp binary_op) {
        return std::accumulate(c.begin(), c.end(), val, binary_op);
    }
    template<class Container>
    inline void remove_duplicate(Container& c) {
        sort(c);
        unique(c);
    }
}
}
