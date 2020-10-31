#include "decl.h"
#include "ILP.h"

bool patchgen::adjust_parameter(PatchParam& param, int variable_index, bool is_increase) {
    // goal: increase/decrease specified variable with minimal change to the rest
    // let: v[i]:=actual variables sought, u[i]:=older vaiables, w[i]:=abs(v[i]-u[i])
    // objective: minimize sum { w[i] }
    int num_sides = param.get_num_sides();

    Eigen::MatrixXd A = get_constraint_matrix(num_sides, param.pattern_id);
    Eigen::VectorXd b = get_constraint_rhs(num_sides, param.pattern_id, param.get_l_permuted());
    int num_variables = A.cols();
    int num_variables_2 = num_variables * 2;

    auto get_variable = [&](PatchParam& param, int index) -> int& {
        return patchgen::get_variable(num_sides, param.pattern_id, param, index);
    };
        
    // auxiliary variables to handle objective of this form:  minimize { abs(x) + abs(y) + ... }
    // actual variables:    v[0], ..., v[n-1]                        n: num_variables
    // auxiliary variables: w[0], ..., w[n-1]
    // appended: v[0], ..., v[n-1], w[0], ..., w[n-1]
    A.conservativeResize(num_sides, num_variables_2);
    A.rightCols(num_variables).setZero();
    ILP ilp(num_variables_2);

    // basic constraints
    ilp.add_constraint(A, EQ, b);

    // constrain the target variable to be increased/decreased
    ilp.add_constraint(Eigen::VectorXd::Unit(num_variables_2, variable_index), is_increase ? GE : LE, get_variable(param, variable_index) + (is_increase ? 1 : -1));

    // constraints for auxiliary variables
    for (int i = 0; i < num_variables; ++i) {
        Eigen::VectorXd row_constr = Eigen::VectorXd::Zero(num_variables_2);
        row_constr[num_variables + i] = -1;
        int old_value = get_variable(param, i);

        row_constr[i] =  1;    ilp.add_constraint(row_constr, LE, old_value);      // constraint:   v[i]-u[i]  <= w[i]
        row_constr[i] = -1;    ilp.add_constraint(row_constr, LE, -old_value);      // constraint: -(v[i]-u[i]) <= w[i]
    }

    // objective
    Eigen::VectorXd row_obj(num_variables_2);
    row_obj << Eigen::VectorXd::Zero(num_variables), Eigen::VectorXd::Ones(num_variables);
    ilp.set_objective(row_obj, false);

    if (!ilp.solve()) return false;

    auto variables = ilp.get_variables();
    for (int i = 0; i < num_variables; ++i)
        get_variable(param, i) = variables[i];

    return true;
}
