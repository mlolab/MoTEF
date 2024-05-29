syms c1 c2 c3 c4 c5 c6;
syms alpha rho n ell tau C;
syms clambda cgamma ceta;



zero = sym(0);
one = sym(1);


clambda = sym(1)/sym(200);
cgamma = sym(1)/sym(200);
ceta = sym(1)/sym(100000);
c = [sym(0.0020); sym(0.000065); sym(0.005); sym(0.0000025); sym(0.01); sym(0.000005)];


zero = sym(0);
one = sym(1);
S = [(one)/(alpha * rho^3*tau*n*ell), zero, zero, zero, zero, zero;
    zero, (one)/(n*ell), zero, zero, zero, zero;
    zero, zero, (ell)/(rho^3*n*tau), zero, zero, zero;
    zero, zero, zero, (one)/(rho*n*ell), zero, zero;
    zero, zero, zero, zero, (ell)/(rho^3*n*tau), zero;
    zero, zero, zero, zero, zero, (one)/(rho*n*ell)];

cc = S * c

lambda = (clambda*alpha^2*rho^6*tau^2)/n;
gamma = cgamma * alpha * rho;
eta = (ceta * alpha * rho^3*tau)/ell;


M = param_matrix_vr(lambda, gamma, eta, n, ell, alpha, rho, C);
b2 = b2_vec_vr(gamma, n, alpha, rho);
qprime = qprime_vec( eta, n, ell);
Mcc = M .* cc.';

MccQprime = [Mcc, qprime];
disp("Perform simplications on the system of linear inequalites by removing the greatest common divisor of each row of MccQprime. Note that the gcd's are positive.")
MSimp = row_rm_gcd(MccQprime)

xx = [one; one; one; one; one; one; -one];

surplus = MSimp * xx
disp("One can check that it sufficies to set alpha, rho, tau, C to be 1, 1, 1, 4 respectively.")
surplus_subs = subs(surplus, [alpha, rho, tau, n, C], [1 1 1 1 4])
