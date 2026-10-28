%% Project 2 - Convection-diffusion equation solver
% -D*Laplacian(u) + beta.grad(u) = 0   in Omega = [0,1]x[0,1]
% u = 1 su Gamma1, u = 0 su Gamma2  
clear; clc; close all;

meshdir = 'mesh2';          
betax = 1; betay = 3;        % velocity field beta = (1,3)^T 
beta = [betax; betay];
%% --- 0. Loading datas and starting computations
% --- Loading mesh and boundary conditions ---
[Nelem, Nodes, triang, coord, NDir, DirNod, DirVal, NeuNod] = input_data2(meshdir);

% --- Computation local basis' coefficents and area ---
[Bloc, Cloc, Area] = localBasis2(Nelem, triang, coord);

% --- Characteristic dimension of meshes ---

h = max(2*sqrt(Area));

% --- RHS: null forcing (f = 0) ---
rhs0 = rhsBuild2(Nodes, coord, Area, triang);

%% --- 0.5. Preliminary check: beta = 0 must reduce to the Poisson problem ---
% Here we only verify that the system solves without
% errors with D=1 and beta=0; the quantitative comparison must be done by also
% running the Project 1 solver on the same mesh with f=0 and the same BC.
A_test = stiffBuild2(Nelem, Nodes, triang, Bloc, Cloc, Area, 1);
B_test = transportBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, 0, 0);
[sysMat_test, rhs_test] = imposeBC2(NDir, DirNod, DirVal, A_test + B_test, rhs0, 1e-10);
u_test = sysMat_test \ rhs_test; 


%% --- 1.: Solution without stabilization for any D ---
D_list = [1, 1e-1, 1e-3]; % Vary the diffusion coefficient D, not delta
sol_noStab = cell(size(D_list));

% Compute the Stiffness matrix with D=1 
A_D1 = stiffBuild2(Nelem, Nodes, triang, Bloc, Cloc, Area, 1); 
B = transportBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, betax, betay);

figure('Name', ' Convection-diffusion without stabilization');
for k = 1:length(D_list)
    D = D_list(k);
    
    % No stabilization (S = 0), let just scale A_D1 for the current D
    sysMat = D * A_D1 + B; 
    
    [sysMat_bc, rhs_bc] = imposeBC2(NDir, DirNod, DirVal, sysMat, rhs0, 1e-10);
    uh = sysMat_bc \ rhs_bc;
    
    sol_noStab{k} = uh; % Save in the correct array for the final table
    
    subplot(2, 2, k);
    trisurf(triang, coord(:,1), coord(:,2), uh, 'EdgeColor', 'k', 'LineWidth', 0.1);
    colormap(gca, parula);
    title(sprintf('D = %.1e, No Stabilizzazione', D));   
    xlabel('x'); ylabel('y'); zlabel('u');
    axis tight; grid on;
    view(120, 25);
end

% extra figure to visualize the solution in the most critic case (D smallest)
subplot(2, 2, 4);
trisurf(triang, coord(:,1), coord(:,2), uh, 'EdgeColor', 'k', 'LineWidth', 0.1);
colormap(gca, parula);
title(sprintf('D = %.1e (alternative view)', D_list(end)));  
xlabel('x'); ylabel('y'); zlabel('u');
axis tight; grid on;
view(35, 25);

%% --- 2. Search for the critical value of D and "pseudo" Peclet number ---
% A robust indicator of numerical oscillation: the exact solution of the
% continuous problem satisfies the maximum principle (0<=u<=1, since the
% boundary data are 0 and 1 and there is no forcing term). If the non-stabilized
% FEM solution exits [0,1] (significantly), we are in
% the presence of numerical oscillations
tol_osc = 1e-2;
D_scan = logspace(0, -3, 40); % from 1 a 1e-3
D_crit = NaN;
for D = D_scan
    A = D * A_D1;
    sysMat = A + B;
    [sysMat_bc, rhs_bc] = imposeBC2(NDir, DirNod, DirVal, sysMat, rhs0, 1e-10);
    uh = sysMat_bc \ rhs_bc;
    if max(uh) > 1 + tol_osc || min(uh) < -tol_osc
        D_crit = D;
        break
    end
end

if ~isnan(D_crit)
    Pe_h = norm(beta) * h / (2 * D_crit);
    fprintf('\nValue of D below which oscillations appear: D_crit = %.4e\n', D_crit);
    fprintf('Characteristic mesh size: h = %.4e\n', h);
    fprintf('Pseudo Peclet number for the meshes: Pe_h = |beta|*h/(2*D_crit) = %.4f\n', Pe_h);
else
    fprintf('\nNo oscillation detected in the explored range of D [%.1e, %.1e].\n', D_scan(end), D_scan(1));
end


%% --- 3. Solution with stabilization Streamline-Upwind Diffusion ---
D_stab = 1e-3;               % critical value for which oscillations occur without stabilization
delta_list = [5e-2, 1e-2, 1e-3];
sol_stab = cell(size(delta_list));

% Base matrix already multiplied by the critical D
A_stabD = D_stab * A_D1; 

figure('Name', ' Convection-diffusion with stabilization SUD');
for k = 1:length(delta_list)
    delta = delta_list(k);
    
    % Construction of the stabilization matrix and system assembly
    S = stabBuild2(Nelem, Nodes, triang, Bloc, Cloc, Area, betax, betay, delta);
    sysMat = A_stabD + B + S;
    
    [sysMat_bc, rhs_bc] = imposeBC2(NDir, DirNod, DirVal, sysMat, rhs0, 1e-10);
    uh = sysMat_bc \ rhs_bc;
    
    sol_stab{k} = uh; 
    
    subplot(2, 2, k);
    trisurf(triang, coord(:,1), coord(:,2), uh, 'EdgeColor', 'k', 'LineWidth', 0.1);
    colormap(gca, parula);
   
    title(sprintf('D = %.1e, delta = %.1e', D_stab, delta));   
    xlabel('x'); ylabel('y'); zlabel('u');
    axis tight; grid on;
view(120, 25);
end

% extra figure: alternative view of the last delta
subplot(2, 2, 4);
trisurf(triang, coord(:,1), coord(:,2), uh, 'EdgeColor', 'k', 'LineWidth', 0.1);
colormap(gca, parula);
title(sprintf('D = %.1e, delta = %.1e (alternative view  alt.)', D_stab, delta_list(end)));  
xlabel('x'); ylabel('y'); zlabel('u');
axis tight; grid on;
view(35, 25);
%% --- Recap Table ---
fprintf('\n%-10s %-12s %-10s %-10s\n', 'D', 'delta', 'max(u)', 'min(u)');
for k = 1:length(D_list)
    fprintf('%-10.2e %-12s %-10.4f %-10.4f\n', D_list(k), '0 (no stab)', max(sol_noStab{k}), min(sol_noStab{k}));
end
for k = 1:length(delta_list)
    fprintf('%-10.2e %-12.2e %-10.4f %-10.4f\n', D_stab, delta_list(k), max(sol_stab{k}), min(sol_stab{k}));
end
