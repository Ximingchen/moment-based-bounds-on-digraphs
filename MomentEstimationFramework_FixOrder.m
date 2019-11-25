function [rho_low, rho_upp, rho_low_refine, rho_upp_refine, w_max_upper] = MomentEstimationFramework_FixOrder(A)
    %% Obtaining the lower bound using moment framework
    n = size(A, 1);
    trA2 = trace(A^2)/n;    % the moments
    trA3 = trace(A^3)/n;
    trA4 = trace(A^4)/n;
    trA5 = trace(A^5)/n;
    
    specRad = max(abs(eig(A)));
    rhol_upper = specRad + 1;
    rhol_lower = 0;
    tol = 1e-3;
    
    Status = 'N';
    while (rhol_upper - rhol_lower > tol) || Status == 'N'
        rho = (rhol_upper + rhol_lower)/2;

        cvx_begin quiet sdp 
            variable u20
            variable u02
            variable u30
            variable u12
            variable u22
            variable u40
            variable u04
            variable u14
            variable u32
            variable u50

            M = [1    0   0 u20 0 u02;...
                 0   u20  0 u30 0 u12; ...
                 0    0  u02 0 u12 0; ...
                 u20 u30  0 u40 0 u22;...
                 0    0  u12 0 u22 0; ...
                 u02 u12  0 u22 0 u04];

            Mx = [0   u20   0   u30   0  u12; ...
                  u20 u30   0   u40   0  u22; ...
                  0    0   u12   0   u22  0 ; ...
                  u30 u40   0   u50   0  u32; ...
                  0    0   u22   0   u32  0 ; ...
                  u12 u22   0   u32   0  u14];

            My = [0    0   u02   0   u12  0; ...
                  0    0   u12   0   u22  0; ...
                  u02 u12   0   u22   0  u04;...
                  0    0   u22   0   u32  0; ...
                  u12 u22   0   u32   0  u14;...
                  0    0   u04   0   u14  0;];
            minimize 0
            subject to
                M >= 0;
                rho*M - Mx >= 0;
                Mx + rho*M >= 0;
                rho*M - My >= 0;
                My + rho*M >= 0;
                u20 >= 0;
                u02 >= 0;
                u22 >= 0;
                u40 >= 0;
                u04 >= 0; 
                u20 - u02 == trA2;
                u30 - 3*u12 == trA3;
                u40 - 6*u22 + u04 == trA4;
                u50 - 10*u32 +5*u14 == trA5;
        cvx_end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhol_upper = rho;
        else
            Status = 'N';
            rhol_lower = rho;
        end
    end
    %disp(['Estimated Lower Bound: ', num2str(rho)]);
    rho_low = rho;

    %% Finding upper bound
    rhou_upper = 20*specRad;
    rhou_lower = specRad - 1;
    tol = 1e-3;
    Status = 'N';
    f = n/(n-1);
    while (rhou_upper - rhou_lower > tol) || Status == 'N'
        rho_upperbound = (rhou_upper + rhou_lower)/2;
        cvx_begin quiet sdp 
            variable u20
            variable u02
            variable u30
            variable u12
            variable u22
            variable u40
            variable u04
            variable u14
            variable u32
            variable u50

            M = [1    -rho_upperbound/(n-1)             0 (n*u20 - rho_upperbound^2)/(n-1)    0         f*u02;...
                 -rho_upperbound/(n-1) (n*u20 - rho_upperbound^2)/(n-1)   0 (n*u30 - rho_upperbound^3)/(n-1)    0         f*u12; ...
                 0    0                              f*u02               0                  f*u12         0; ...
                 (n*u20 - rho_upperbound^2)/(n-1) (n*u30 - rho_upperbound^3)/(n-1)  0 (n*u40 - rho_upperbound^4)/(n-1) 0 f*u22;...
                 0    0  f*u12 0 f*u22 0; ...
                 f*u02 f*u12  0 f*u22 0 f*u04];

            Mx = [-rho_upperbound/(n-1)   (n*u20 - rho_upperbound^2)/(n-1)   0   (n*u30 - rho_upperbound^3)/(n-1)   0  f*u12; ...
                  (n*u20 - rho_upperbound^2)/(n-1) (n*u30 - rho_upperbound^3)/(n-1)   0   (n*u40 - rho_upperbound^4)/(n-1)   0  f*u22; ...
                  0    0   f*u12   0   f*u22  0 ; ...
                  (n*u30 - rho_upperbound^3)/(n-1) (n*u40 - rho_upperbound^4)/(n-1)   0  (n*u50 - rho_upperbound^5)/(n-1)   0  f*u32; ...
                  0    0   f*u22   0   f*u32  0 ; ...   
                  f*u12 f*u22   0   f*u32   0  f*u14];

            My = f*[0    0   u02   0   u12  0; ...
                  0    0   u12   0   u22  0; ...
                  u02 u12   0   u22   0  u04;...
                  0    0   u22   0   u32  0; ...
                  u12 u22   0   u32   0  u14;...
                  0    0   u04   0   u14  0;];

            minimize 0
            subject to
                M >= 0;
                rho_upperbound*M - Mx >= 0;
                Mx + rho_upperbound*M >= 0;
                rho_upperbound*M - My >= 0;
                My + rho_upperbound*M >= 0;
                %rho_upperbound*M - My >= 0;
                %My + rho_upperbound*M >= 0;
                u20 >= 0;
                u02 >= 0;
                u22 >= 0;
                u40 >= 0;
                u04 >= 0; 
                u20 - u02 == trA2;
                u30 - 3*u12 == trA3;
                u40 - 6*u22 + u04 == trA4;
                u50 - 10*u32 +5*u14 == trA5;
        cvx_end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhou_lower = rho_upperbound;
        else
            Status = 'N';
            rhou_upper = rho_upperbound;
        end
    end
    disp(['Estimated Upper Bound: ', num2str(rho_upperbound)]);
    rho_upp = rho_upperbound;
    
    %% Upper and lower bound for Imaginary parts
    B = (A - A')/(2*1i);
    r_up = max(abs(eig(B)));
    r_lim = 8*r_up; 
    tol = 1e-5;
    
    Traces_hat = zeros(5, 2);
    Traces_hat(1, 1) = trace(B^0);
    Traces_hat(3, 1) = trace(B^2);
    Traces_hat(5, 1) = trace(B^4);
    Traces_hat(2, 2) = Traces_hat(5, 1);
    Traces_hat(4, 2) = trace(B^6);
    while (abs(r_up - r_lim) > tol)
        r_avg = (r_up + r_lim)/2;
        [H0_up, H1_up] = buildHankelUpperBound(B, 3, r_avg, 1, Traces_hat);
        [~,check0] = chol(H0_up);
        [~,check1] = chol(r_avg*H0_up - H1_up);
        [~,check2] = chol(r_avg*H0_up + H1_up);
        checker = (check0  == 0 & check1 == 0 & check2 == 0)
        if checker
            r_up = r_avg; % both matrices are positive semidefinite
        else % infeasible
            r_lim = r_avg;
        end
    end
    w_max_upper = r_lim;
    %% Refinements - Lower Bound
	rhol_upper = specRad + 1;
    rhol_lower = 0;
    Status = 'N';
    while (rhol_upper - rhol_lower > tol) || Status == 'N'
        rho = (rhol_upper + rhol_lower)/2;

        cvx_begin quiet sdp 
            variable u20
            variable u02
            variable u30
            variable u12
            variable u22
            variable u40
            variable u04
            variable u14
            variable u32
            variable u50

            M = [1    0   0 u20 0 u02;...
                 0   u20  0 u30 0 u12; ...
                 0    0  u02 0 u12 0; ...
                 u20 u30  0 u40 0 u22;...
                 0    0  u12 0 u22 0; ...
                 u02 u12  0 u22 0 u04];

            Mx = [0   u20   0   u30   0  u12; ...
                  u20 u30   0   u40   0  u22; ...
                  0    0   u12   0   u22  0 ; ...
                  u30 u40   0   u50   0  u32; ...
                  0    0   u22   0   u32  0 ; ...
                  u12 u22   0   u32   0  u14];

            My = [0    0   u02   0   u12  0; ...
                  0    0   u12   0   u22  0; ...
                  u02 u12   0   u22   0  u04;...
                  0    0   u22   0   u32  0; ...
                  u12 u22   0   u32   0  u14;...
                  0    0   u04   0   u14  0;];
            minimize 0
            subject to
                M >= 0;
                rho*M - Mx >= 0;
                Mx + rho*M >= 0;
                w_max_upper*M - My >= 0;
                My + w_max_upper*M >= 0;
                
                rho*M - My >= 0;
                My + rho*M >= 0;
                u20 >= 0;
                u02 >= 0;
                u22 >= 0;
                u40 >= 0;
                u04 >= 0; 
                u20 - u02 == trA2;
                u30 - 3*u12 == trA3;
                u40 - 6*u22 + u04 == trA4;
                u50 - 15*u32 +5*u14 == trA5;
        cvx_end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhol_upper = rho;
        else
            Status = 'N';
            rhol_lower = rho;
        end
    end
    %disp(['Estimated Lower Bound: ', num2str(rho)]);
    rho_low_refine = rho;

    %% Refinement - Upper Bound 
    rhou_upper = 300;
    rhou_lower = specRad - 1;
    Status = 'N';
    f = n/(n-1);
    while (rhou_upper - rhou_lower > tol) || Status == 'N'
        rho_upperbound = (rhou_upper + rhou_lower)/2;
        cvx_begin quiet sdp 
            variable u20
            variable u02
            variable u30
            variable u12
            variable u22
            variable u40
            variable u04
            variable u14
            variable u32
            variable u50

            M = [1    -rho_upperbound/(n-1)                                 0 (n*u20 - rho_upperbound^2)/(n-1)    0         f*u02;...
                 -rho_upperbound/(n-1)   (n*u20 - rho_upperbound^2)/(n-1)   0 (n*u30 - rho_upperbound^3)/(n-1)    0         f*u12; ...
                 0    0                              f*u02               0                  f*u12         0; ...
                 (n*u20 - rho_upperbound^2)/(n-1) (n*u30 - rho_upperbound^3)/(n-1)  0 (n*u40 - rho_upperbound^4)/(n-1) 0 f*u22;...
                 0    0  f*u12 0 f*u22 0; ...
                 f*u02 f*u12  0 f*u22 0 f*u04];

            Mx = [-rho_upperbound/(n-1)   (n*u20 - rho_upperbound^2)/(n-1)   0   (n*u30 - rho_upperbound^3)/(n-1)   0  f*u12; ...
                  (n*u20 - rho_upperbound^2)/(n-1) (n*u30 - rho_upperbound^3)/(n-1)   0   (n*u40 - rho_upperbound^4)/(n-1)   0  f*u22; ...
                  0    0   f*u12   0   f*u22  0 ; ...
                  (n*u30 - rho_upperbound^3)/(n-1) (n*u40 - rho_upperbound^4)/(n-1)   0  (n*u50 - rho_upperbound^5)/(n-1)   0  f*u32; ...
                  0    0   f*u22   0   f*u32  0 ; ...   
                  f*u12 f*u22   0   f*u32   0  f*u14];

            My = f*[0    0   u02   0   u12  0; ...
                  0    0   u12   0   u22  0; ...
                  u02 u12   0   u22   0  u04;...
                  0    0   u22   0   u32  0; ...
                  u12 u22   0   u32   0  u14;...
                  0    0   u04   0   u14  0;];

            minimize 0
            subject to
                M >= 0;
                rho_upperbound*M - Mx >= 0;
                Mx + rho_upperbound*M >= 0;
                w_max_upper*M - My >= 0;
                My + w_max_upper*M >= 0;

                %rho_upperbound*M - My >= 0;
                %My + rho_upperbound*M >= 0;
                u20 >= 0;
                u02 >= 0;
                u22 >= 0;
                u40 >= 0;
                u04 >= 0; 
                u20 - u02 == trA2;
                u30 - 3*u12 == trA3;
                u40 - 6*u22 + u04 == trA4;
                u50 - 15*u32 +5*u14 == trA5;
        cvx_end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhou_lower = rho_upperbound;
        else
            Status = 'N';
            rhou_upper = rho_upperbound;
        end
    end
    %disp(['Estimated Upper Bound: ', num2str(rho_upperbound)]);
    rho_upp_refine = rho_upperbound;
end