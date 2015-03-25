%referenced by crosslink.m
function dc = cross_fun(t,c,n,d, kon, koff, s,g, c0)
    M = 2*pi;
    x = M*d*(n/g/(n+M) + (n-1)/(M+n-1));
    c1 = c(1);
    c2 = c(2);
    dc = [2*kon*(c0-c1-c2)*(s*2/g-c1-2*c2) - koff*c1 + koff*c2 - kon*c1*(s*2/g-c1-2*c2)*x;
        - koff*c2 + kon*c1*(s*2/g-c1-2*c2)*x];
end