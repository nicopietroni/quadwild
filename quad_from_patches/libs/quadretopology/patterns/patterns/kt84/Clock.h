#pragma once

#include <ctime>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>

namespace kt84 {

class ClkStart;
class ClkData
{
    friend ClkStart;
    std::map<std::string, clock_t> m_ticks;
    std::string m_clkName;
    bool m_showPercent;
public:
    ClkData(const std::string& clkName = "", bool showPercent = true)
        : m_clkName(clkName)
        , m_showPercent(showPercent)
    {}
    void setClkName(const std::string& clkName) {
        m_clkName = clkName;
    }
    std::string getClkName() const {
        return m_clkName;
    }
    void clearTicks() { m_ticks.clear(); }
    clock_t getTick(const std::string& name) const {
        auto p = m_ticks.find(name);
        return p == m_ticks.end() ? 0 : p->second;
    }
    void setShowPercent(bool showPercent) {
        m_showPercent = showPercent;
    }
    void print() {
        double total = 0.0;
        for (auto i = m_ticks.begin(); i != m_ticks.end(); ++i)
            total += static_cast<double>(i->second);
        std::cout << "[";
        if (!m_clkName.empty())
            std::cout << m_clkName << " | ";
        for (auto i = m_ticks.begin(); i != m_ticks.end(); ++i) {
            std::cout << i->first << ":";
            if (m_showPercent) {
                double percent = total == 0.0 ? 0.0 : (100 * i->second / total);
                std::cout << std::setprecision(3) << percent << "%, ";
            } else {
                std::cout << (1000 * i->second / static_cast<double>(CLOCKS_PER_SEC)) << "msec, ";
            }
        }
        std::cout << "total:" << (1000 * total / CLOCKS_PER_SEC) << "msec]\n";
    }
};

class ClkStart
{
    ClkData* m_clkData;
    std::string m_currentTick;
    bool m_autoPrint;
    clock_t m_clk;
public:
    ClkStart(ClkData* clkData, const std::string& currentTick = "", bool autoPrint = false, bool clearTicks = false)
        : m_clkData(clkData)
        , m_currentTick(currentTick)
        , m_autoPrint(autoPrint)
        , m_clk(clock())
    {
        if (m_autoPrint)
            std::cout << "clock started: " << m_clkData->getClkName() << std::endl;
        if (clearTicks)
            m_clkData->clearTicks();
    }
    void setCurrentTick(const std::string& tick) {
        clock_t clk_old = m_clk;
        m_clk = clock();
        if (m_currentTick != "")
            m_clkData->m_ticks[m_currentTick] += m_clk - clk_old;
        m_currentTick = tick;
    }
    void print() {
        clock_t clk_old = m_clk;
        m_clk = clock();
        m_clkData->m_ticks[m_currentTick] += m_clk - clk_old;
        m_clkData->print();
    }
    ~ClkStart() {
        clock_t clk_old = m_clk;
        m_clk = clock();
        m_clkData->m_ticks[m_currentTick] += m_clk - clk_old;
        if (m_autoPrint) {
            m_clkData->print();
            std::cout << "clock stopped: " << m_clkData->getClkName() << std::endl;
        }
    }
};
class ClkSimple
{
    std::string m_message;
    clock_t m_clk;
public:
    ClkSimple(const std::string& message = "")
        : m_message(message)
        , m_clk(clock())
    {
        std::cout << "Start " << m_message << "...";
    }
    ~ClkSimple()
    {
        std::cout << "done: " << (1000 * (clock() - m_clk) / CLOCKS_PER_SEC) << "msec\n";
    }
    void print() {
        std::cout << "current: " << (1000 * (clock() - m_clk) / CLOCKS_PER_SEC) << "msec\n";
    }
};

class ClkLoop {
    int m_tick;
    int m_tickMax;
    int m_printIntervalMilliSec;
    std::string m_name;
    clock_t m_clk, m_clkInit;
public:
    ClkLoop()
        : m_tick(0)
        , m_tickMax(0)
        , m_printIntervalMilliSec(0)
        , m_clk(0)
        , m_clkInit(0)
    {}
    ClkLoop(int tickMax, std::string name = "", int printIntervalMilliSec = 1000) {
        init(tickMax, name, printIntervalMilliSec);
    }
    void init(int tickMax, std::string name = "", int printIntervalMilliSec = 1000) {
        m_tick = 0;
        m_tickMax = tickMax;
        m_printIntervalMilliSec = printIntervalMilliSec;
        m_name = name;
        m_clk = m_clkInit = clock();
    }
    void tick() {
        ++m_tick;
        clock_t clk_new = clock();
        if (m_tick == m_tickMax || m_printIntervalMilliSec < (clk_new - m_clk) * 1000 / CLOCKS_PER_SEC) {
            double percent = (100. * m_tick) / m_tickMax;
            int step = static_cast<int>(percent / 10);
            std::cout << m_name << " |";
            for (int i = 0; i < 10; ++i)
                std::cout << (i < step ? "#" : "-");
            std::cout << "| " << percent << "% ";
            if (m_tick < m_tickMax)
                std::cout << "     \r";
            else
                std::cout << ((clk_new - m_clkInit) / CLOCKS_PER_SEC) <<"sec\n";
            m_clk = clk_new;
        }
    }
    void done() {     // i'm done!
        if (m_tick == m_tickMax)
            return;
        m_tick = m_tickMax - 1;
        tick();
    }
};

}

