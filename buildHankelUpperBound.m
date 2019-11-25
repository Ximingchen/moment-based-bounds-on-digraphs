function [H0, H1] = buildHankelUpperBound(A, order, rho, provided, Traces)
    %if trace(A^2)/sum(sum(A)) ~= 1
    %    error('Graph is directed, cannot use this function');
    %end
    if nargin <= 3
        provided = 0;
        Traces = [];
    end
	n = size(A, 1);
    c0 = zeros(order+1,1); r0 = c0;
    c1 = zeros(order+1,1); r1 = c1;
    if provided
        for i = 1:order + 1
            c0(i) = (Traces(i, 1) - ((-1)^(i-1) + 1)* rho^(i-1))/(n-2);
            r0(i) = (Traces(i, 2)- ((-1)^(i-1+order) + 1)*rho^(i-1 + order))/(n-2);
            c1(i) = (Traces(i+1, 1) - ((-1)^(i) + 1)*rho^i)/(n-2);
            r1(i) = (Traces(i+1, 2) - ((-1)^(i+order) + 1)*rho^(i+order))/(n-2);
        end
    else
        for i = 1:order + 1
            c0(i) = (trace(A^(i-1)) - ((-1)^(i-1) + 1)* rho^(i-1))/(n-2);
            r0(i) = (trace(A^(i-1+order))- ((-1)^(i-1+order) + 1)*rho^(i-1 + order))/(n-2);
            c1(i) = (trace(A^(i)) - ((-1)^(i) + 1)*rho^i)/(n-2);
            r1(i) = (trace(A^(i+order)) - ((-1)^(i+order) + 1)*rho^(i+order))/(n-2);
        end
    end
    H0 = hankel(c0, r0);
    H1 = hankel(c1, r1);
end