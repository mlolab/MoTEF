syms c1 c2 c3 c4 c5 c6;
syms alpha rho n L tau C;
syms clambda cgamma ceta;
zero = sym(0);
one = sym(1);


clambda = sym(1)/sym(200);
cgamma = sym(1)/sym(200);
ceta = sym(1)/sym(100000);



c = [sym(0.0020); sym(0.000065); sym(0.05); sym(0.0000025); sym(0.09); sym(0.000005)];



S = [(one)/(n^2*L), zero, zero, zero, zero, zero;
    zero, ( tau)/(n*L), zero, zero, zero, zero;
    zero, zero, (L)/(rho^3 * n*tau), zero, zero, zero;
    zero, zero, zero, (tau)/(rho*n*L), zero, zero;
    zero, zero, zero, zero, (L)/(rho^3*n*tau), zero;
    zero, zero, zero, zero, zero, (tau)/(rho*n*L)];

cc = S * c

lambda = (clambda*alpha*tau*rho^3);
gamma = cgamma * alpha * rho;
eta = (ceta * alpha *tau *rho^3)/L;


M = param_matrix(lambda, gamma, eta, n, L, alpha, rho, C);
b2 = b2_vec(gamma, n, alpha, rho);
qprime = qprime_vec( eta, n, L);
Mcc = M .* cc.';


MccQprime = [Mcc, qprime];

disp("Perform simplications on the system of linear inequalites by removing the greatest common divisor of each row of MccQprime. Note that the gcd's are positive.")
MSimp = row_rm_gcd(MccQprime)

xx = [one; one; one; one; one; one; -one];

surplus = MSimp * xx
disp("One can check that it sufficies to set alpha, rho, tau, C to be 1, 1, 1, 4 respectively.")
surplus_subs = subs(surplus, [alpha, rho, tau, C], [1 1 1 4])
