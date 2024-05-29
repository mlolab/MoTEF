function M = param_matrix_pl(lambda, gamma, eta, n, L, mu, alpha, rho, C)
    zero = sym(0);
    M = [(lambda - mu* eta), zero, zero, zero, zero, zero;
         zero, (lambda - mu* eta), zero, (- (sym(6) * lambda^2) / (alpha)), zero, (- (sym(12) * lambda^2)/ (gamma * rho));
        %  (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), ((alpha)/sym(2) - (sym(6)*gamma^2*C)/alpha), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha),(- (sym(6) * gamma * C)/rho), ( - (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
         (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), ((alpha)/sym(4) - mu* eta), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha),(- (sym(6) * gamma * C)/rho), ( - (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
        %  zero, zero, zero, ((alpha)/sym(2) - (sym(6)*gamma^2*C)/alpha), zero, (- (sym(6) * gamma * C)/rho);
         zero, zero, zero, ((alpha)/sym(4)  - mu* eta), zero, (- (sym(6) * gamma * C)/rho);
         (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), (- (sym(6)*gamma^2*C)/alpha), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha), ((gamma * rho)/sym(2)  - mu* eta), (- (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
         (-(sym(3) * L^2 * n * eta^2)/ lambda), (-(sym(3) * L^2 * eta^2)/ lambda), (-(sym(6) * eta^2)/alpha), (-(sym(6) * gamma^2 * C)/alpha - (sym(36)*lambda^2*eta^2*L^2)/alpha), (- (sym(6) * eta^2)/(gamma *rho)), ((gamma * rho)/ sym(2) -(sym(36)*eta^2*lambda^2*L^2)/(gamma*rho)  - mu* eta);
         (-(sym(3) * L^2 * n^2 * eta^2)/ lambda), (-(sym(3) * L^2 * n * eta^2)/ lambda), (-(sym(6) * eta^2 * n)/alpha),(- (sym(36)*lambda^2*eta^2*L^2)/alpha), zero, (-(sym(36) * eta^2 * gamma*L^2)/rho)
         ];
end