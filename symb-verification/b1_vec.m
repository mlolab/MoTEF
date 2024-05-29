function b= b1_vec(lambda, gamma, eta, n, L, alpha, rho)
    b = [
        (sym(3) * L^2 * n^2 * eta^2) / lambda;
        (sym(3) * L^2 * n * eta^2)/lambda;
        (sym(6) * eta^2 * n) / alpha;
        (sym(36) * eta^2 * lambda^2 * L^2) / alpha; 
        sym(0);
        (sym(36) * eta^2 * gamma * L^2) / rho
    ];
end