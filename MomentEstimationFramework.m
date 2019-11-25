function output = MomentEstimationFramework(A, r, shape, optimize, disp_info)
    %% input checking
    switch nargin
        case 1
            r = 5;
            shape = 'square';
            optimize = 0;
            disp_info = false;            
        case 2
            shape = 'square';
            optimize = 0;
            disp_info = false;            
        case 3
            optimize = 0;
            disp_info = false;
        case 4
            disp_info = false;
    end
    %% Parameter settings
    n = size(A, 1);
    d = 2;
    scaled_traces = zeros(r, 1);
    for i = 1:r
        scaled_traces(i) = trace(A^i)/n;    % the moments
    end  
    specRad = max(abs(eig(A)));   
    rho_low = computeLower(specRad, scaled_traces, r, d, shape); %, w);
    rho_upp = computeUpper(specRad, scaled_traces, n, r, d, shape); %, shape);
    symmetrized_order = 3;
    w = computeOmegaMax(A, symmetrized_order);
    rho_sym = computeRhoSym(A, symmetrized_order);
    if optimize
        rho_low_refined = computeLower(specRad, scaled_traces, r, d, shape, w, disp_info);
        rho_upp_refined = computeUpper(specRad, scaled_traces, n, r, d, shape, w, disp_info);
    else
        rho_low_refined = rho_low;
        rho_upp_refined = rho_upp;
    end
    output.rho_low = rho_low;
    output.rho_upp = rho_upp;
    output.w = w;
    output.rho_low_refined = rho_low_refined;
    output.rho_upp_refined = rho_upp_refined;
    output.rho_symmetrized_upper = rho_sym;
end

function rho_low = computeLower(specRad, scaled_traces, r, d, shape, w, disp_info)
    switch nargin
        case 4
            shape = 'square';
            use_omega = false;
            w = 100;
            disp_info = false;
        case 5
            use_omega = false;
            w = 100;
            disp_info = false;
        case 6
            use_omega = true;
            disp_info = false;
        case 7
            use_omega = true;
        otherwise
            use_omega = false;
    end
    
    r2 = floor(r/2);                    % the monomial degree 
    r_bino = nchoosek(d + r, r);        % the number of moments considered in the optimization
    r2_bino = nchoosek(d + r2, r2);     % the dimension of moment matrix
    indices = generate_multiindex(r);   % this sets the order of moments considered to be at most r
    assert(size(indices, 1) == r_bino, 'Error: dimension mismatch between multisequence index and C(n+r)(r)');
    rhol_upper = specRad;           % bisection parameters for lower bound
    rhol_lower = 0;                     
    tol = 1e-3;                         % tolerence to stop the bisection process
    tol_precision = 1/(10^3);           % precision to deal with numerical issues in equality constraints in optimization
    Status = 'N';                       % default status of the optimization program
    
    r_v = floor(r/2) - 1;               % the monomial degree for the circle case
    r_v_bino = nchoosek(d + r_v, r_v);  % the dimension of the localizing matrix in the circle case
    if mod(r, 2) == 0   % when we have truncated sequence up to an even order, we have to remove one dimension to construct localizing matrix
        r_b = r2 - 1;
        r_b_bino = nchoosek(d+r_b, r_b);
    else    % otherwise since the degree of linear function is 1, then we floor(r/2)*2 + 1 = r
        r_b = r2;
        r_b_bino = r2_bino;
    end
    if  disp_info
        disp('Calculating lower bound...');
        disp(['Using ', shape, ' support']);
    end
    while (rhol_upper - rhol_lower > tol) 
        rho = (rhol_upper + rhol_lower)/2;
        if disp_info
            disp(['upper: ',num2str(rhol_upper), ' lower: ',num2str(rhol_lower)]);
        end
        if shape == 'square'
            cvx_begin quiet sdp 
                variable u(r_bino)
                M = cvx(zeros(r2_bino, r2_bino));
                Mx = cvx(zeros(r_b_bino, r_b_bino));
                My = cvx(zeros(r_b_bino, r_b_bino));
                y = cvx(zeros(r, 1));
                % set the trace sums
                for i = 1:r
                    for s = 0:floor(i/2)
                        [~, idx] = ismember(indices,[i - 2*s, 2*s], 'rows');
                        idx = find(idx);
                        y(i) = y(i) + nchoosek(i, 2*s)*(-1)^s*u(idx);
                    end
                end
                % set the Moment matrix
                for i = 1:r2_bino
                    for j = i:r2_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta), 'rows');
                        idx = find(idx);
                        M(i,j) = u(idx);
                    end
                end
                M = triu(M, 1) + triu(M)';
                % set the localizing matrices
                Mb = M(1:r_b_bino, 1:r_b_bino); 
                % an auxillary matrix to consider when r is an even number
                % it is used to match the dimension with Mx and My.
                for i = 1:r_b_bino
                    for j = i:r_b_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta + [1, 0]), 'rows');
                        idx = find(idx);
                        Mx(i,j) = u(idx);
                    end
                end
                Mx = triu(Mx, 1) + triu(Mx)';
                % set My
                for i = 1:r_b_bino
                    for j = i:r_b_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta + [0, 1]), 'rows');
                        idx = find(idx);
                        My(i,j) = u(idx);
                    end
                end
                My = triu(My, 1) + triu(My)';
                minimize 0
                subject to
                    % constraints due to feasibility of moments
                    M >= 0;
                    rho*Mb - Mx >= 0;
                    Mx + rho*Mb >= 0;
                    if use_omega
                        w*Mb - My >= 0;
                        My + w*Mb >= 0;
                    else
                        rho*Mb - My >= 0;
                        My + rho*Mb >= 0;    
                    end
                    
                    % spectrum constraints
                    for i = 1:r_bino
                        alpha = indices(i, :);
                        if mod(alpha(2), 2) == 1    % when b is an odd number, the mixed moment is zero
                            u(i) >= 0;
                            u(i) <= tol_precision;
                        else
                            if mod(alpha(1), 2) == 0 && mod(alpha(2), 2) == 0
                                if (alpha(1)+alpha(2) > 0)
                                    u(i) >= 0;
                                else
                                    u(i) == 1;
                                end
                            end
                        end
                    end
                    for i = 1:r
                       y(i) - scaled_traces(i) <= tol_precision;
                       y(i) - scaled_traces(i) >= -tol_precision;
                    end
            cvx_end
            
        else   % the case when the support is a circle
            if mod(r, 2) == 0 %even number
                r_box = floor(r/2) - 1;
            else
                r_box = floor(r/2);
            end
            r_box_bino = nchoosek(d + r_box, r_box);
            cvx_begin quiet sdp 
                variable u(r_bino)
                M = cvx(zeros(r2_bino, r2_bino));
                M2 = cvx(zeros(r_v_bino, r_v_bino));
                Mx2 = cvx(zeros(r_v_bino, r_v_bino));
                My2 = cvx(zeros(r_v_bino, r_v_bino));
                My_box = cvx(zeros(r_box_bino,r_box_bino));
                y = cvx(zeros(r, 1));
                % set the trace sums
                for i = 1:r
                    for s = 0:floor(i/2)
                        [~, idx] = ismember(indices,[i - 2*s, 2*s], 'rows');
                        idx = find(idx);
                        y(i) = y(i) + nchoosek(i, 2*s)*(-1)^s*u(idx);
                    end
                end
                % set the Moment matrix
                for i = 1:r2_bino
                    for j = i:r2_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta), 'rows'); idx = find(idx);
                        M(i,j) = u(idx);
                    end
                end
                M = triu(M, 1) + triu(M)';
                M_box = M(1:r_box_bino, 1:r_box_bino);
                % set the localizing matrices
                for i = 1:r_v_bino
                    for j = i:r_v_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta), 'rows'); idx = find(idx);
                        M2(i,j) = u(idx);
                    end
                end
                M2 = triu(M2, 1) + triu(M2)';
                for i = 1:r_v_bino
                    for j = i:r_v_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta + [2, 0]), 'rows'); idx = find(idx);
                        Mx2(i,j) = u(idx);
                    end
                end
                Mx2 = triu(Mx2, 1) + triu(Mx2)';
                % set My
                for i = 1:r_v_bino
                    for j = i:r_v_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta + [0, 2]), 'rows'); idx = find(idx);
                        My2(i,j) = u(idx);
                    end
                end
                My2 = triu(My2, 1) + triu(My2)';
                for i = 1:r_box_bino
                    for j = i:r_box_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices,(alpha + beta + [0, 1]), 'rows');
                        idx = find(idx);
                        My_box(i,j) = u(idx);
                    end
                end
                My_box = triu(My_box, 1) + triu(My_box)';
                
                minimize 0
                subject to
                    % constraints due to feasibility of moments
                    M >= 0;
                    rho^2*M2 - Mx2 - My2 >= 0;
                    if use_omega
                        w*M_box - My_box >= 0;
                        My_box + w*M_box >= 0;
                    end
                    % spectrum constraints
                    for i = 1:r_bino
                        alpha = indices(i, :);
                        if mod(alpha(2), 2) == 1
                            u(i) == 0;
                        else
                            if mod(alpha(1), 2) == 0 && mod(alpha(2), 2) == 0
                                if (alpha(1)+alpha(2) > 0)
                                    u(i) >= 0;
                                else
                                    u(i) == 1;
                                end
                            end
                        end
                    end
                    for i = 1:r
                        y(i) - scaled_traces(i) <= tol_precision;
                        y(i) - scaled_traces(i) >= -tol_precision;
                    end
                cvx_end
        end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhol_upper = rho;
        else
            Status = 'N';
            rhol_lower = rho;
        end
    end
    rho_low = rho;
    if disp_info
        disp(['Estimated Lower Bound: ', num2str(rho_low)]);
    end
end

function rho_upp = computeUpper(specRad, scaled_traces, n, r, d, shape, w, disp_info) % shape) 
    switch nargin
        case 4
            d = 2;
            shape = 'square';
            use_omega = false;
            w = 100;
            disp_info = false;
        case 5
            shape = 'square';
            use_omega = false;
            w = 1;
            disp_info = false;
        case 6
            use_omega = false;
            disp_info = false;
        case 7
            use_omega = true;
            disp_info = false;
        otherwise
            if nargin == 8
                use_omega = true;
            else
                use_omega = false;
            end
            if nargin <= 3
                disp_info = false;
            end
    end
    
	rhou_upper = 10*specRad;    % bisection parameter for upper bound
    rhou_lower = specRad;   % bisection parameter
    factor = n/(n-1);           % the scaling factor since the largest atom is removed
    r2 = floor(r/2);            
    r_bino = nchoosek(d + r, r);
    r2_bino = nchoosek(d + r2, r2);
    indices = generate_multiindex(r);   % this sets the order of moments considered to be at most r
    assert(size(indices, 1) == r_bino, 'Error: dimension mismatch between multisequence index and C(n+r)(r)');
    tol = 1e-3;                 % tolerence to stop the bisection process
    Status = 'N';               % default status
    tol_precision = 1/(2*n);   % tolerence to avoid numerical issues in SDP
    
	r_v = floor(r/2) - 1;       % the monomial degree for the case of circle
    r_v_bino = nchoosek(d + r_v, r_v);  % the dimension of the localizing matrix in the case of a circle support
    
    if mod(r, 2) == 0
        r_b = r2 - 1;
        r_b_bino = nchoosek(d + r_b, r_b);
    else
        r_b = r2;
        r_b_bino = r2_bino;
    end
    
    if  disp_info
        disp('Computing upper bound ...');
        disp(['Using ', shape, ' support']);
    end
    while (rhou_upper - rhou_lower > tol) 
        rho_upperbound = (rhou_upper + rhou_lower)/2;
        if disp_info
            disp(['upper: ',num2str(rhou_upper), ' lower: ',num2str(rhou_lower)]);
        end
        if shape == 'square'
            cvx_begin quiet sdp 
                variable u(r_bino)
                variable tmp_obj
                M = cvx(zeros(r2_bino, r2_bino));
                Mx = cvx(zeros(r_b_bino, r_b_bino));
                My = cvx(zeros(r_b_bino, r_b_bino));
                y = cvx(zeros(r, 1));
                %% set the trace sums
                for i = 1:r
                    for s = 0:floor(i/2)
                        [~, idx] = ismember(indices,[i - 2*s, 2*s], 'rows');
                        idx = find(idx);
                        y(i) = y(i) + nchoosek(i, 2*s)*(-1)^s*u(idx);
                    end
                end
                % set the Moment matrix
                for i = 1:r2_bino
                    for j = i:r2_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        moment = alpha + beta;
                        [~, idx] = ismember(indices, moment, 'rows');
                        idx = find(idx);
                        if moment(2) ~= 0
                            M(i,j) = factor*u(idx);
                        else
                            M(i,j) = (n*u(idx) - rho_upperbound^moment(1))/(n-1);
                        end
                    end
                end
                M = triu(M, 1) + triu(M)';
                Mb = M(1:r_b_bino, 1:r_b_bino);
                
                % set the localizing matrices
                for i = 1:r_b_bino
                    for j = i:r_b_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        moment = alpha + beta + [1, 0];
                        [~, idx] = ismember(indices, moment, 'rows');
                        idx = find(idx);
                        if moment(2) ~= 0
                            Mx(i,j) = factor*u(idx);
                        else
                            Mx(i,j) = (n*u(idx) - rho_upperbound^moment(1))/(n-1);
                        end
                    end
                end
                Mx = triu(Mx, 1) + triu(Mx)';
                % set My
                for i = 1:r_b_bino
                    for j = i:r_b_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices, (alpha + beta + [0, 1]), 'rows');
                        idx = find(idx);
                        My(i,j) = factor.*u(idx);
                    end
                end
                My = triu(My, 1) + triu(My)';

                minimize 0
                subject to
                    % constraints due to feasibility of moments
                    M >= 0;
                    rho_upperbound*Mb - Mx >= 0;
                    Mx + rho_upperbound*Mb >= 0;
                    if use_omega
                        w*Mb - My >= 0;
                        My + w*Mb >= 0;
                    else
                        rho_upperbound*Mb - My >= 0;
                        My + rho_upperbound*Mb >= 0;
                    end
                    % spectrum constraints
                    for i = 1:r_bino
                        alpha = indices(i, :);
                        if mod(alpha(2), 2) == 1
                            u(i) == 0;
                        else
                            if mod(alpha(1), 2) == 0 && mod(alpha(2), 2) == 0
                                if (alpha(1)+alpha(2) > 0)
                                    u(i) >= 0;
                                else
                                    u(i) == 1;
                                end
                            end
                        end
                    end
                    for i = 1:r
                       y(i) - scaled_traces(i) <= tol_precision;
                       y(i) - scaled_traces(i) >= -tol_precision;
                    end
            cvx_end
        else          
            if mod(r, 2) == 0 %even number
                r_box = floor(r/2) - 1;
            else
                r_box = floor(r/2);
            end
            r_box_bino = nchoosek(d + r_box, r_box);
            cvx_begin quiet sdp 
                variable u(r_bino)
                variable tmp_obj
                M = cvx(zeros(r2_bino, r2_bino));
                Mx2 = cvx(zeros(r_v_bino, r_v_bino));
                My2 = cvx(zeros(r_v_bino, r_v_bino));
                %Mx = cvx(zeros(r2_bino, r2_bino));
                My_box = cvx(zeros(r_box_bino, r_box_bino));
                y = cvx(zeros(r, 1));
                %% set the trace sums
                for i = 1:r
                    for s = 0:floor(i/2)
                        [~, idx] = ismember(indices,[i - 2*s, 2*s], 'rows');
                        idx = find(idx);
                        y(i) = y(i) + nchoosek(i, 2*s)*(-1)^s*u(idx);
                    end
                end
                % set the Moment matrix
                for i = 1:r2_bino
                    for j = i:r2_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        moment = alpha + beta;
                        [~, idx] = ismember(indices, moment, 'rows');
                        idx = find(idx);
                        if moment(2) ~= 0
                            M(i,j) = factor*u(idx);
                        else
                            M(i,j) = (n*u(idx) - rho_upperbound^moment(1))/(n-1);
                        end
                    end
                end
                M = triu(M, 1) + triu(M)';
                M_box = M(1:r_box_bino,1:r_box_bino);
                % set the localizing matrices
                M2 = M(1:r_v_bino, 1:r_v_bino);
                % set Mx2
                for i = 1:r_v_bino
                    for j = i:r_v_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        moment = (alpha + beta + [2, 0]);
                        [~, idx] = ismember(indices, moment, 'rows'); idx = find(idx);
                        if moment(2) ~= 0
                            Mx2(i,j) = factor*u(idx);
                        else
                            Mx2(i,j) = (n*u(idx) - rho_upperbound^moment(1))/(n-1);
                        end
                    end
                end
                Mx2 = triu(Mx2, 1) + triu(Mx2)';
                % set My2
                for i = 1:r_v_bino
                    for j = i:r_v_bino
                        alpha = indices(i, :); beta = indices(j, :);
                        moment = (alpha + beta + [0, 2]);
                        [~, idx] = ismember(indices, moment, 'rows');
                        idx = find(idx);
                        My2(i,j) = factor.*u(idx);
                    end
                end
                My2 = triu(My2, 1) + triu(My2)';
                for i = 1:r_box_bino
                    for j = i:r_box_bino
                        alpha = indices(i, :);
                        beta = indices(j, :);
                        [~, idx] = ismember(indices, (alpha + beta + [0, 1]), 'rows');
                        idx = find(idx);
                        My_box(i,j) = factor.*u(idx);
                    end
                end
                My_box = triu(My_box, 1) + triu(My_box)';
                minimize 0
                subject to
                    % constraints due to feasibility of moments
                    M >= 0;
                    rho_upperbound^2*M2 - Mx2 - My2 >= 0;
                    if use_omega
                        w*M_box - My_box >= 0;
                        My_box + w*M_box >= 0;
                    end
                    % spectrum constraints
                    for i = 1:r_bino
                        alpha = indices(i, :);
                        if mod(alpha(2), 2) == 1
                            u(i) == 0;
                        else
                            if mod(alpha(1), 2) == 0 && mod(alpha(2), 2) == 0
                                if (alpha(1)+alpha(2) > 0)
                                    u(i) >= 0;
                                else
                                    u(i) - 1 <= tol_precision;
                                    u(i) - 1 >= -tol_precision;
                                end
                            end
                        end
                    end
                    for i = 1:r
                       y(i) - scaled_traces(i) <= tol_precision;
                       y(i) - scaled_traces(i) >= -tol_precision;
                    end
            cvx_end            
        end
        if strcmp(cvx_status, 'Solved')
            Status = 'Y';
            rhou_lower = rho_upperbound;
        else
            Status = 'N';
            rhou_upper = rho_upperbound;
        end
    end
    if disp_info
        disp(['Estimated Upper Bound: ', num2str(rho_upperbound)]);
    end
    rho_upp = rho_upperbound;
end
