function A = GenerateGraphs(n, w, Mode, isDigraph, isSimple)
    if nargin ~= 5
        isDigraph = 1;
        isSimple = 1;
    end
    if nargin <3 
        Mode = 'ER';
    end
    
    switch Mode
        case 'ER'
            assert(length(w) == 1, 'Should be a single probability');
            A = rand(n, n) < w;
        case 'CL'
            assert(size(w, 1) == n, 'Size of weight vector should be consistent with the dimension');
            w_in = w(:, 1);
            w_out = w(:, 2);
            
            NormFactor = sum(w_in);
            assert(sum(w_in) == sum(w_out), 'In-degree must equal to out-degree');
            
            P = w_in * w_out';
            A = rand(n, n) < (P/NormFactor);
        % case 'SBM'
        otherwise
            error('No such option');
    end
    
    if ~isDigraph
        A = triu(A);
        A = A + A';
    end
    if isSimple
        A = A - diag(diag(A));
    end
end