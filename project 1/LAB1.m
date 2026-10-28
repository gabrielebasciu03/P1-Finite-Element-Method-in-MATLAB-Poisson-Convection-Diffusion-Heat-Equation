%% Project 1 - Finite element code for 2D piecewise linear Galerkin

% -div(grad u)(x) = f(x) with Dirichlet and Neumann BCs
clear; clc;
meshes = {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4'};
errors = zeros(length(meshes), 1);
h = zeros(length(meshes), 1);
tol   = 1e-8;
maxit = 1000;
allResHist = cell(length(meshes),1);
allRealResHist = cell(length(meshes),1);
for m = 1:length(meshes)
    meshdir = meshes{m};
    
    % Loading datas
    [Nelem, Nodes, triang, coord, NDir, DirNod, DirVal, ureal, force,NeuNod] = input_data(meshdir);
    
    %% --- 1. Computation local basis' coefficents and area ---
    [Bloc, Cloc, Area] = localBasis(Nelem, triang, coord);
    
    %% --- 2. Construction stiffness matrix ---
    stiffMat = stiffBuild(Nelem, Nodes, triang, Bloc, Cloc, Area);
   
    %% --- 3. Construction load vector (RHS) and imposing Neumann ---
    rhs = rhsBuild(Nodes, coord, Area, triang, force,NeuNod);
     sysMat = stiffMat;

    %% 4. Imposing Dirichlet conditions (Penalization Method) ---
    [sysMat, rhs] = imposeBC(NDir, DirNod, DirVal, sysMat, rhs);

    %% --- 5. Solution linear system and definition h 
    uh = sysMat \ rhs;
    [uh_pcg, resHist, realResHist] = pcgSolve(sysMat, rhs, tol, maxit);
    fprintf('%s: |uh - uh_pcg| = %.3e, iterations PCG = %d\n', ...
    meshes{m}, norm(uh - uh_pcg), length(resHist));

    allResHist{m}     = resHist;
    allRealResHist{m} = realResHist;
    h(m) = max(2*sqrt(Area));  

    %% --- 6. Computation L2 error ---
    ur = ureal(coord(:,1), coord(:,2));
    nodal_area = calculateNodalArea(Nodes, triang, Area);
    errors(m) = sqrt(sum(((uh - ur).^2) .* nodal_area));
    
   
    
    %% --- 7. Plots ---
    if m == length(meshes)
        figure(1);
        trisurf(triang, coord(:,1), coord(:,2), uh);
        title('Numerical Solution uh');
        figure(2);
        trisurf(triang, coord(:,1), coord(:,2), ur);
        title('Real Solution u');

        figure(3);
        
        trisurf(triang, coord(:,1), coord(:,2), ur - uh);
        title('Error Function u - uh');
    end
end
fprintf('%-8s %-12s %-12s %-10s\n', 'Mesh', 'h', ' L2 error', 'ratio');
for m = 1:length(meshes)
    if m == 1
        ratio = NaN;
    else
        ratio = errors(m-1)/errors(m);
    end
    fprintf('%-8s %-12.4e %-12.4e %-10.4f\n', meshes{m}, h(m), errors(m), ratio);
end


%  Creation figure
figure(4);

% Plot of the error (blue line with circles)
loglog(h, errors, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
       'DisplayName', 'Errore L_2 calcolato (\epsilon_h)');
hold on;
%  Creation of the reference's theoretical line O(h^2)
% To make the comparison clearer let's multiply the line for a constant C
% so that the reference line is starting from the same starting point 
% of the error.
C = errors(1) / (h(1)^2);
reference_line = C * (h.^2);

% Plot of the referene pendence (blach line)
loglog(h, reference_line, '--k', 'LineWidth', 1.5, ...
       'DisplayName', 'Theoretical pendence O(h^2)');


grid on;

grid minor; 

xlabel('mesh''s characteristic dimension (h)');
ylabel(' L_2 error norm (\epsilon_h)');
title('Convergence of the FEM error (P_1)');
legend('Location', 'northwest');
hold off;

% --- Convergence graphic PCG for all the meshes (semi-log) ---
figure(5);
hold on;
colors = lines(length(meshes));
for m = 1:length(meshes)
    plot(1:length(allResHist{m}), allResHist{m}, '-', ...
        'Color', colors(m,:), 'LineWidth', 1.5, ...
        'DisplayName', meshes{m});
end
grid on;
xlabel('iteration PCG (k)');
ylabel('log_{10}(||r_k|| / ||b||)');
title('PCG convergence for all the meshes (preconditioned IC(0))');
legend('Location', 'northeast');
hold off;

% --- Comparison ''interior'' residual PCG vs real residual b - A*x_k ---
% Comparison only on two meshes (mesh 0 and mesh 4)
% per non affollare troppo il grafico.
figure(6);
subplot(1,2,1);
plot(1:length(allResHist{1}), allResHist{1}, '-o', 'DisplayName', ' PCG residual (interior)');
hold on;
plot(1:length(allRealResHist{1}), allRealResHist{1}, '--x', 'DisplayName', 'real residual b-Ax_k');
grid on; xlabel('iteration k'); ylabel('log_{10}(relative residual)');
title([' residual comparison - ', meshes{1}]);
legend('Location','northeast');
hold off;

subplot(1,2,2);
plot(1:length(allResHist{end}), allResHist{end}, '-o', 'DisplayName', 'PCG residual (interior)');
hold on;
plot(1:length(allRealResHist{end}), allRealResHist{end}, '--x', 'DisplayName', 'real residual b-Ax_k');
grid on; xlabel('Iteration k'); ylabel('log_{10}(relative residual)');
title([' residual comparison - ', meshes{end}]);
legend('Location','northeast');
hold off;
p_global = polyfit(log(h), log(errors), 1);
fprintf('\nGlobal convergence order (least-squares fit): %.4f  (expected ~2 for P1 FEM)\n', p_global(1));
