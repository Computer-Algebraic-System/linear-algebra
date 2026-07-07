#pragma once

inline std::map<algebra::Variable, double> tensor::solve_linear_system(const std::vector<algebra::Equation>& equations, const Method method) {
    std::map<algebra::Variable, double> res;
    std::vector<double> values;

    for (const algebra::Equation& equation : equations) {
        GLOBAL_FORMATTING << equation << std::endl;
    }
    const auto [A, X, B] = Matrix<algebra::Variable>::from_equations(equations);

    switch (method) {
    case Method::GAUSS:
        values = A.augment(B).gauss_elimination();
        break;

    case Method::CRAMER:
        values = A.cramer_rule(B);
        break;
    }
    const int size = values.size();

    for (int i = 0; i < size; i++) {
        res.emplace(X[i, 0], values[i]);
        GLOBAL_FORMATTING << algebra::Equation(X[i, 0], values[i]) << std::endl;
    }
    return res;
}
