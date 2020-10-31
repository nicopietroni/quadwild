#pragma once
#include <vector>
#include <list>
#include <set>

namespace kt84 {
    template<typename T_to, typename T_from>
	std::vector<T_to> container_cast(const std::vector<T_from>& container) {
        std::vector<T_to> result;
        result.reserve(container.size());
        for (auto& e : container)
            result.push_back(T_to(e));
        return result;
    }
    template<typename T_to, typename T_from>
	std::list<T_to> container_cast(const std::list<T_from>& container) {
        std::list<T_to> result;
        for (auto& e : container)
            result.push_back(T_to(e));
        return result;
    }
    template<typename T_to, typename T_from>
	std::set<T_to> container_cast(const std::set<T_from>& container) {
        std::set<T_to> result;
        for (auto& e : container)
            result.insert(T_to(e));
        return result;
    }
}
