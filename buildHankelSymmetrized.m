function [H0, H1] = buildHankelSymmetrized(A, order, rho, Traces)
    n = size(A, 1);
	c0 = zeros(order+1,1); r0 = c0;
    c1 = zeros(order+1,1); r1 = c1;
    for i = 1:order + 1
        c0(i) = (Traces(i, 1) -  rho^(i-1))/(n-1);
        r0(i) = (Traces(i, 2)-   rho^(i-1 + order))/(n-1);
        c1(i) = (Traces(i+1, 1) - rho^i)/(n-1);
        r1(i) = (Traces(i+1, 2) - rho^(i+order))/(n-1);
        
    end
    H0 = hankel(c0, r0);
    H1 = hankel(c1, r1);
end