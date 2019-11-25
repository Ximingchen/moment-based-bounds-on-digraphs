% beta is the power-law exponent, i.e., the number of nodes with degree k
% is proportional to k^-beta
% Delta is the maximum expected degree
% d is the average expected degree
function w = PowerLawCoef(beta, Delta, d, n)
    assert(beta >= 3, 'beta should be larger or equal to 3');
    i0 = n*(d*(beta - 2)/ (Delta*(beta - 1)))^(beta - 1);
    c = (beta - 2)/(beta - 1) * d*n^(1/(beta - 1));

    w = zeros(n, 1);
    for i = 1:n
        w(i) = c*(i0 + i)^(-1/(beta - 1));
    end
end