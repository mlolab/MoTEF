function MS = row_rm_gcd(M)
    [Mrows, Mcols] = size(M);  % Get the size of M
    MS = sym(zeros(size(M)));  % Initialize MS to have same size as M
    for i = 1:Mrows
        row_gcd = multi_gcd(M(i,:));   % Calculate GCD
        row_gcd
        MS(i,:) = M(i,:) / row_gcd;     % Divide row by GCD
    end
end