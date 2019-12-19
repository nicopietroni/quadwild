#pragma once
#include <lp_lib.h>
//#include "gurobi_c++.h"
#include <Eigen/Core>

struct ILP {
    int num_variables;
    lprec* ptr;
    
    ILP(int num_variables)
        : num_variables(num_variables)
        , ptr(::make_lp(0, num_variables))
    {
        ::set_verbose(ptr, SEVERE);
        // set all variables to be integer
        for (int i = 0; i < num_variables; ++i)
            if (!::set_int(ptr, i + 1, TRUE)) assert(false);
        // default ojective
        set_objective(Eigen::VectorXd::Ones(num_variables), true);
    }
    ~ILP() { ::delete_lp(ptr); }
    void add_constraint(const Eigen::VectorXd& row, int constr_type, double rhs) const {
        assert(row.size() == num_variables);
        Eigen::VectorXd tmp;
        tmp.resize(num_variables + 1);
        tmp[0] = 0;
        tmp.tail(num_variables) = row;
        if (!::add_constraint(ptr, &tmp[0], constr_type, rhs)) assert(false);
    }
    void add_constraint(const Eigen::MatrixXd& rows, int constr_type, const Eigen::VectorXd& rhs) const {
        for (int i = 0; i < rows.rows(); ++i)
            add_constraint(rows.row(i), constr_type, rhs[i]);
    }
    void set_objective(const Eigen::VectorXd& row, bool is_maxim) const {
        assert(row.size() == num_variables);
        Eigen::VectorXd tmp;
        tmp.resize(num_variables + 1);
        tmp[0] = 0;
        tmp.tail(num_variables) = row;
        if (!::set_obj_fn(ptr, &tmp[0])) assert(false);
        (is_maxim ? ::set_maxim : ::set_minim)(ptr);
    }
    bool solve() const { return ::solve(ptr) <= SUBOPTIMAL; }
    void refresh() {            // workaround for weird error when re-solving
        auto ptr_old = ptr;
        ptr = ::copy_lp(ptr);
        ::delete_lp(ptr_old);
    }
    Eigen::VectorXi get_variables() const {
        Eigen::VectorXd variables = Eigen::VectorXd::Zero(num_variables);
        if (!::get_variables(ptr, &variables[0])) assert(false);
        return variables.cast<int>();
    }
};

