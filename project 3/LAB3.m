%% Project 3 - nparabolic_fem
%  Finite Element Method (P1 Galerkin) + theta-method in time
%  for the 2D parabolic problem
%
%    du/dt = Delta u + f(x,y,t)     in Omega = [-1,1]^2
%    u = 0                          on d(Omega), for all t in [0,T]
%    u(x,y,0) = u0(x,y)
%
%  Manufactured solution:
%     u(x,y,t) = exp(-t) * sin(pi*x) * cos(pi/2*y)
%     f(x,y,t) = (5/4*pi^2 - 1) * exp(-t) * sin(pi*x) * cos(pi/2*y)
%
%  theta-method:
%     (1/dt * P + theta*A) u^{n+1} =
%         (1/dt * P - (1-theta)*A) u^n + theta*b^{n+1} + (1-theta)*b^n
%
%  theta = 0    -> explicit Euler   (conditionally stable)
%  theta = 0.5  -> Crank-Nicolson   (unconditionally stable, O(dt^2))
%  theta = 1    -> implicit Euler   (unconditionally stable, O(dt)).

clear all; close all; clc;

meshes  = {'mesh0','mesh1','mesh2','mesh3','mesh4'};

%% ---- manufactured solution / forcing --------------------------------
uex   = @(x,y,t) exp(-t) .* sin(pi.*x) .* cos(pi/2.*y);
force = @(x,y,t) (5*pi^2/4 - 1) .* exp(-t) .* sin(pi.*x) .* cos(pi/2.*y);

T = 1;

%% ---  1. Spatial convergence (theta = 1 and theta = 0.5)---
%  Fixed, very small dt (so the temporal error is negligible) while h
%  is refined mesh0 -> mesh4. Expected slope ~ 2 for both theta.

fprintf('\n');
fprintf('Spatial convergence Explicit Euler and Crank Nicolson');
fprintf('\n');

dt_small   = 1e-4;                 % small enough to isolate space error
thetaListA = [1, 0.5];
hvals      = zeros(numel(meshes),1);
errA       = zeros(numel(meshes), numel(thetaListA));

for im = 1:numel(meshes)
    meshdir = meshes{im};
    for it = 1:numel(thetaListA)
        theta = thetaListA(it);
        [errST, ~, ~, hloc] = solveParabolic(meshdir, theta, dt_small, T, uex, force);
        hvals(im)        = hloc;
        errA(im, it)      = errST;
        fprintf('mesh=%-6s theta=%-4.2f  h=%.4e  err=%.4e\n', ...
                 meshdir, theta, hloc, errST);
    end
end

figure('Name','Spatial convergence');
loglog(hvals, errA(:,1), 'bo-', 'LineWidth',1.5, 'MarkerSize',7); hold on;
loglog(hvals, errA(:,2), 'rs-', 'LineWidth',1.5, 'MarkerSize',7);
co1 = polyfit(log(hvals), log(errA(:,1)), 1);
co2 = polyfit(log(hvals), log(errA(:,2)), 1);
loglog(hvals, exp(polyval(co1, log(hvals))), 'b--');
loglog(hvals, exp(polyval(co2, log(hvals))), 'r--');
legend(sprintf('\\theta=1, slope=%.2f', co1(1)), ...
       sprintf('\\theta=0.5, slope=%.2f', co2(1)), ...
       'fit \theta=1','fit \theta=0.5','Location','northwest');
xlabel('h'); ylabel('L^2(0,T;L^2) error'); grid on;
title('Spatial convergence (dt fixed, small)');
fprintf('Spatial order (theta=1)  : %.4f  (expected ~2)\n', co1(1));
fprintf('Spatial order (theta=0.5): %.4f  (expected ~2)\n', co2(1));

%% --- 2. Temporal convergence (theta = 1 and theta = 0.5)---
%  Fixed fine mesh, dt refined. Expected slope ~1 for theta=1,
%  slope ~2 for theta=0.5 (Crank-Nicolson).

fprintf('\n');
fprintf('Temporal convergence');
fprintf('\n');

meshFine   = 'mesh4';                       % fine enough: space error << time error
dtList     = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625];
errB       = zeros(numel(dtList), numel(thetaListA));

for id = 1:numel(dtList)
    dt = dtList(id);
    for it = 1:numel(thetaListA)
        theta = thetaListA(it);
        [errST, ~, ~, ~] = solveParabolic(meshFine, theta, dt, T, uex, force);
        errB(id, it) = errST;
        fprintf('theta=%-4.2f  dt=%.5f  err=%.4e\n', theta, dt, errST);
    end
end

figure('Name','Temporal convergence');
loglog(dtList, errB(:,1), 'bo-', 'LineWidth',1.5, 'MarkerSize',7); hold on;
loglog(dtList, errB(:,2), 'rs-', 'LineWidth',1.5, 'MarkerSize',7);
cot1 = polyfit(log(dtList), log(errB(:,1)'), 1);
cot2 = polyfit(log(dtList), log(errB(:,2)'), 1);
loglog(dtList, exp(polyval(cot1, log(dtList))), 'b--');
loglog(dtList, exp(polyval(cot2, log(dtList))), 'r--');
legend(sprintf('\\theta=1, slope=%.2f', cot1(1)), ...
       sprintf('\\theta=0.5, slope=%.2f', cot2(1)), ...
       'fit \theta=1','fit \theta=0.5','Location','northwest');
xlabel('\Delta t'); ylabel('L^2(0,T;L^2) error'); grid on;
title(sprintf('Temporal convergence (mesh=%s fixed)', meshFine));
fprintf('Temporal order (theta=1)  : %.4f  (expected ~1)\n', cot1(1));
fprintf('Temporal order (theta=0.5): %.4f  (expected ~2)\n', cot2(1));
%% ---2.5.  Self-convergence (Richardson) for temporal order ---

%  Compares u_h^{dt} vs u_h^{dt/2} on the SAME mesh: the spatial error
%  cancels out, isolating the temporal error and its true order,
%  even when it falls below the spatial error floor.

fprintf('\n');
fprintf('Self-convergence for time');
fprintf('\n');

meshFine   = 'mesh4';                     % same fine mesh as Part B
dtListSC   = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125];

errSC = zeros(numel(dtListSC), numel(thetaListA));

for it = 1:numel(thetaListA)
    theta = thetaListA(it);
    for id = 1:numel(dtListSC)
        dt = dtListSC(id);

        [~, ~, uh_dt,   ~, nodeArea] = solveParabolic(meshFine, theta, dt,   T, uex, force);
        [~, ~, uh_dt2,  ~, ~       ] = solveParabolic(meshFine, theta, dt/2, T, uex, force);

        diffSol = uh_dt - uh_dt2;
        errSC(id, it) = sqrt(dot(diffSol.^2, nodeArea));

        fprintf('theta=%-4.2f  dt=%.6f  ||u^dt - u^{dt/2}||=%.4e\n', ...
                 theta, dt, errSC(id, it));
    end
end

figure('Name','Self-convergence (Richardson) in time');
loglog(dtListSC, errSC(:,1), 'bo-', 'LineWidth',1.5, 'MarkerSize',7); hold on;
loglog(dtListSC, errSC(:,2), 'rs-', 'LineWidth',1.5, 'MarkerSize',7);
cosc1 = polyfit(log(dtListSC), log(errSC(:,1)'), 1);
cosc2 = polyfit(log(dtListSC), log(errSC(:,2)'), 1);
loglog(dtListSC, exp(polyval(cosc1, log(dtListSC))), 'b--');
loglog(dtListSC, exp(polyval(cosc2, log(dtListSC))), 'r--');
legend(sprintf('\\theta=1, slope=%.2f', cosc1(1)), ...
       sprintf('\\theta=0.5, slope=%.2f', cosc2(1)), ...
       'fit \theta=1','fit \theta=0.5','Location','northwest');
xlabel('\Delta t'); ylabel('||u_h^{\Delta t} - u_h^{\Delta t/2}||_{L^2}'); grid on;
title(sprintf('Self-convergence in time (mesh=%s fixed)', meshFine));
fprintf('Self-convergence order (theta=1)  : %.4f  (expected ~1)\n', cosc1(1));
fprintf('Self-convergence order (theta=0.5): %.4f  (expected ~2)\n', cosc2(1));

%% --- 3. theta = 0 (explicit Euler): stability-limited dt ---
%  
%  Estimate dt_max from the generalized eigenvalue problem
%  A x = lambda P x  ->  dt < 2/lambda_max(P^-1 A).
%  Skip meshes for which the number of required steps blows up
%  (as suggested in the assignment, typically mesh3/mesh4).

fprintf('\n');
fprintf('Implicit Euler')
fprintf('\n');

safety      = 0.9;                 % safety factor on the stability bound
maxSteps    = 20000;                % refuse meshes needing more steps than this
hvalsC      = [];
errC        = [];

for im = 1:numel(meshes)
    meshdir = meshes{im};

    [Nelem,Nodes,triang,coord,~,DirNod] = input_data3(meshdir);
    [Bloc, Cloc, Area] = localBasis3(Nelem, triang, coord);
    A = stiffBuild3(Nelem, Nodes, triang, Bloc, Cloc, Area);
    P = massBuild3(Nelem, Nodes, triang, Area);

    lambda_max = eigs(A, P, 1, 'lm');
    dt_stab    = safety * 2 / lambda_max;
    nsteps     = ceil(T/dt_stab);

    h_est = sqrt(2*max(Area));

    if nsteps > maxSteps
        fprintf('mesh=%-6s h=%.4e  dt_stab=%.3e  steps=%d  -> SKIPPED (too many steps)\n', ...
                 meshdir, h_est, dt_stab, nsteps);
        continue;
    end

    dt_used = T/nsteps;   % exact divisor of T, <= dt_stab
    [errST, ~, ~, hloc] = solveParabolic(meshdir, 0, dt_used, T, uex, force);

    fprintf('mesh=%-6s h=%.4e  dt_stab=%.3e  dt_used=%.3e  steps=%d  err=%.4e\n', ...
             meshdir, hloc, dt_stab, dt_used, nsteps, errST);

    hvalsC(end+1)  = hloc; 
    errC(end+1)    = errST; 
end

if numel(hvalsC) >= 2
    coC = polyfit(log(hvalsC), log(errC), 1);
    figure('Name','theta=0 convergence (stability-limited)');
    loglog(hvalsC, errC, 'ko-', 'LineWidth',1.5, 'MarkerSize',7); hold on;
    loglog(hvalsC, exp(polyval(coC, log(hvalsC))), 'k--');
    legend(sprintf('\\theta=0 error, slope=%.2f', coC(1)), 'fit', ...
           'Location','northwest');
    xlabel('h'); ylabel('L^2(0,T;L^2) error'); grid on;
    title('\theta=0 (explicit Euler): convergence at the stability limit');
    fprintf('theta=0 order : %.4f\n', coC(1));
else
    fprintf('Not enough usable meshes for theta=0 to fit a convergence order.\n');
end



