function q=qprime_vec( eta, n, L)
    zero = sym(0);
    q = [
        eta / (n^2);
        zero;
        zero;
        zero;
        (eta * L^2) / n;
        zero;
        % (eta^2 *L)/sym(2) - eta/sym(2)
        -eta/sym(4)
    ];
end