function M = param_matrix(lambda, gamma, eta, n, L, alpha, rho, C)
    zero = sym(0);
    M = [lambda, zero, zero, zero, zero, zero;
         zero, lambda, zero, (- (sym(6) * lambda^2) / (alpha)), zero, (- (sym(12) * lambda^2)/ (gamma * rho));
        %  (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), ((alpha)/sym(2) - (sym(6)*gamma^2*C)/alpha), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha),(- (sym(6) * gamma * C)/rho), ( - (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
         (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), ((alpha)/sym(4)), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha),(- (sym(6) * gamma * C)/rho), ( - (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
        %  zero, zero, zero, ((alpha)/sym(2) - (sym(6)*gamma^2*C)/alpha), zero, (- (sym(6) * gamma * C)/rho);
         zero, zero, zero, ((alpha)/sym(4)), zero, (- (sym(6) * gamma * C)/rho);
         (-(sym(3) * L^2 * n * gamma^2 * C)/ lambda), (-(sym(3) * L^2  * gamma^2 * C)/ lambda), (- (sym(6)*gamma^2*C)/alpha), (-(sym(36) * lambda^2 * gamma^2*L^2*C)/alpha), ((gamma * rho)/sym(2)), (- (sym(36) * gamma * lambda^2 * L^2 * C)/rho);
         (-(sym(3) * L^2 * n * eta^2)/ lambda), (-(sym(3) * L^2 * eta^2)/ lambda), (-(sym(6) * eta^2)/alpha), (-(sym(6) * gamma^2 * C)/alpha - (sym(36)*lambda^2*eta^2*L^2)/alpha), (- (sym(6) * eta^2)/(gamma *rho)), ((gamma * rho)/ sym(2) -(sym(36)*eta^2*lambda^2*L^2)/(gamma*rho));
         (-(sym(3) * L^2 * n^2 * eta^2)/ lambda), (-(sym(3) * L^2 * n * eta^2)/ lambda), (-(sym(6) * eta^2 * n)/alpha),(- (sym(36)*lambda^2*eta^2*n* L^2)/alpha), zero, (-(sym(36) * eta^2 * n* gamma*L^2)/rho)
         ];
end