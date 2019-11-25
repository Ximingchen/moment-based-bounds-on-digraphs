function w_max_upper = computeOmegaMax(A, order)
    B = (A - A')/(2*1i);
    r_up = max(abs(eig(B)));
    r_lim = 8*r_up; 
    tol = 1e-5;
    
    Traces = zeros(order + 2, 2);
    for k = 1:order + 2
        if mod(k - 1, 2) == 0
            Traces(k, 1) = trace(B^(k-1));
        end
        if mod(k-1+order, 2) == 0
            Traces(k, 2) = trace(B^(k-1+order));
        end
    end
    while (abs(r_up - r_lim) > tol)
        r_avg = (r_up + r_lim)/2;
        [H0_up, H1_up] = buildHankelUpperBound(B, order, r_avg, 1, Traces);
        [~,check0] = chol(H0_up);
        [~,check1] = chol(r_avg*H0_up - H1_up);
        [~,check2] = chol(r_avg*H0_up + H1_up);
        if (check0  == 0 && check1 == 0 && check2 == 0)
            r_up = r_avg; % both matrices are positive semidefinite
        else % infeasible
            r_lim = r_avg;
        end
    end
    w_max_upper = r_lim;
end