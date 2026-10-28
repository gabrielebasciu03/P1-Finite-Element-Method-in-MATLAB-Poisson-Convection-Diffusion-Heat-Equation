function [errL2spacetime, tgrid, uh_final, h_est, nodeArea] = ...
          solveParabolic(meshdir, theta, dt, T, uex, force)
% SOLVERPARABOLIC, solve the parabolic problem with the theta-method on a given mesh,
% and return the L2(0,T;L2(Omega)) norm of the error, plus the final
% solution and nodeArea (needed for self-convergence / Richardson checks).

    [Nelem, Nodes, triang, coord, ~, DirNod] = input_data3(meshdir);
    x = coord(:,1);
    y = coord(:,2);

    [Bloc, Cloc, Area] = localBasis3(Nelem, triang, coord);
    A = stiffBuild3(Nelem, Nodes, triang, Bloc, Cloc, Area);
    P = massBuild3(Nelem, Nodes, triang, Area);
    nodeArea = nodeAreaAssemble(Nelem, Nodes, triang, Area);

    h_est = sqrt(2*max(Area));

    N  = max(1, round(T/dt));
    dt = T/N;                       % make dt an exact divisor of T
    tgrid = linspace(0, T, N+1);

    % Initial condition (exact solution at t=0)
    u_n = uex(x, y, 0);

    % Assemble the (constant, since dt is fixed) LHS matrix and impose BC
    LHS = (1/dt)*P + theta*A;
    LHS = imposeBC3(Nodes, DirNod, LHS);
    dLHS = decomposition(LHS, 'lu');   % reuse the factorization at every step

    RHSmat = (1/dt)*P - (1-theta)*A;

    err2_time = zeros(N+1,1);   % squared spatial L2 error at each time level
    err2_time(1) = 0;           % exact at t = 0

    for n = 0:N-1
        tn   = tgrid(n+1);
        tnp1 = tgrid(n+2);

        bn   = P * force(x, y, tn);
        bnp1 = P * force(x, y, tnp1);

        rhs = RHSmat*u_n + theta*bnp1 + (1-theta)*bn;
        rhs(DirNod) = 0;            % homogeneous Dirichlet BC

        u_np1 = dLHS \ rhs;

        uerr = uex(x, y, tnp1) - u_np1;
        err2_time(n+2) = dot(uerr.^2, nodeArea);

        u_n = u_np1;
    end

    errL2spacetime = sqrt(trapz(tgrid, err2_time));
    uh_final = u_n;
end
