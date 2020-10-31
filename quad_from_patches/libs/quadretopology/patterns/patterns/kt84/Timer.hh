#pragma once
#include <chrono>
#include <unordered_map>

namespace kt84 {
    template<class Rep, class Period>
    struct Timer {
        Timer(bool auto_start = false) : is_active() {
            if (auto_start) start();
        }
        void reset() {
            duration = std::chrono::duration<Rep, Period>::zero();
            is_active = false;
        }
        void start() {
            if (is_active) return;
            clk_before = std::chrono::high_resolution_clock::now();
            is_active = true;
        }
        Rep stop() {
            if (is_active) {
                auto clk_after = std::chrono::high_resolution_clock::now();
                duration += clk_after - clk_before;
                is_active = false;
            }
            return count();
        }
        Rep count() const {
            return duration.count();
        }
        template <typename Func>
        Rep start_stop(Func func) {
            start();
            func();
            return stop();
        }
        std::chrono::high_resolution_clock::time_point clk_before;
        std::chrono::duration<Rep, Period> duration;
        bool is_active;
    };
    typedef Timer<float, std::milli> Timerf_milli;
    typedef Timer<float, std::micro> Timerf_micro;
    typedef Timer<float, std::nano > Timerf_nano ;
    typedef Timer<double, std::milli> Timerd_milli;
    typedef Timer<double, std::micro> Timerd_micro;
    typedef Timer<double, std::nano > Timerd_nano ;
    
    template<class Rep, class Period>
    struct NamedTimers : public std::unordered_map<const char*, Timer<Rep, Period>> {};
    typedef NamedTimers<float, std::milli> NamedTimersf_milli;
    typedef NamedTimers<float, std::micro> NamedTimersf_micro;
    typedef NamedTimers<float, std::nano > NamedTimersf_nano ;
    typedef NamedTimers<double, std::milli> NamedTimersd_milli;
    typedef NamedTimers<double, std::micro> NamedTimersd_micro;
    typedef NamedTimers<double, std::nano > NamedTimersd_nano ;
}
