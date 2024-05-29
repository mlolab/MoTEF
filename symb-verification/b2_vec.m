function b= b2_vec(gamma, n, alpha, rho)
    b = [
        n;
        sym(2) * n;
        sym(0);
        (sym(6) * n) /alpha;
        sym(0);
        (sym(6) * n) / (gamma * rho)
    ];
end