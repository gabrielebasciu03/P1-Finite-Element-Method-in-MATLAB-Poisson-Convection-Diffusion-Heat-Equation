 function [x, resHistory, realResHistory] = pcgSolve(A, b, tol, maxit)
% PCGSOLVE  Preconditioned Conjugate Gradient with Incomplete Cholesky (zero fill-in).
% Tracks the relative residual norm at each iteration for convergence plots.
%
% INPUT:
%   A        - sparse SPD system matrix
%   b        - right hand side vector
%   tol      - relative residual tolerance (e.g. 1e-8)
%   maxit    - maximum number of iterations
%
% OUTPUT:
%   x              - computed solution
%   resHistory     - log10(||r_k|| / ||b||) at each iteration k (PCG "internal" residual)
%   realResHistory  - log10(||b - A*x_k|| / ||b||) at each iteration (true residual, for comparison)

    n = size(A,1);
    x = zeros(n,1);

    % --- Incomplete Cholesky preconditioner with zero fill-in ---
    % ichol requires a matrix with positive diagonal; our stiffness matrix
    % after Dirichlet penalty is SPD, so this should apply directly.
    opts.type = 'nofill';
    L = ichol(A, opts);

    r = b - A*x;
    normb = norm(b);
    if normb == 0
        normb = 1; % avoid division by zero
    end

    z = L' \ (L \ r);
    p = z;
    rz_old = r'*z;

    resHistory = zeros(maxit,1);
    realResHistory = zeros(maxit,1);

    for k = 1:maxit
        Ap = A*p;
        alpha = rz_old / (p'*Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;

        relres = norm(r)/normb;
        resHistory(k) = log10(relres);
        realResHistory(k) = log10(norm(b - A*x)/normb);

        if relres < tol
            resHistory = resHistory(1:k);
            realResHistory = realResHistory(1:k);
            break;
        end

        z = L' \ (L \ r);
        rz_new = r'*z;
        beta = rz_new / rz_old;
        p = z + beta*p;
        rz_old = rz_new;
    end
end
