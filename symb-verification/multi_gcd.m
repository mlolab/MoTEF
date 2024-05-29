function g = multi_gcd(vec)
    g = vec(1);
    for i = 2:length(vec)
        g = gcd(g, vec(i));
    end
end