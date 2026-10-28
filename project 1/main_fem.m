% ==========================================================================
% Finite element code for 2D piecewise linear Galerkin solution of
%   -div(grad u) = f   in Omega = [-1,1]x[-1,1]
%   u = 0                          on y = +-1   (Dirichlet, GammaD)
%   grad(u).nu = qN                on x = +-1   (Neumann,   GammaN)
%
% Numerical Methods for Differential Equations - Lab Project Part I
% ==========================================================================
clear; close all; clc;

meshdirs = {'mesh0','mesh1','mesh2','mesh3','mesh4'};
nlev = numel(meshdirs);

err   = zeros(nlev,1);
hval  = zeros(nlev,1);
nodes_v = zeros(nlev,1);

% store results at an intermediate refinement level for plotting
plot_level = 3;   % mesh2  (1-based index -> meshdirs{3} = 'mesh2')
plot_data = struct();

for lev = 1:nlev
    meshdir = meshdirs{lev};

    % ---- read mesh & problem data --------------------------------------
    [Nelem,Nodes,triang,coord,NDir,DirNod,DirVal,NNeu,NeuNod,ureal,force] = ...
        input_data(meshdir);

    % ---- local basis-function coefficients (Algorithm 1) ---------------
    [Bloc,Cloc,Area] = localBasis(Nelem,triang,coord);

    % ---- assemble stiffness matrix (Algorithm 2) ------------------------
    Stiff = stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area);

    % ---- assemble load vector (domain + Neumann terms) -------------------
    h = 2/(sqrt(Nodes)-1); % uniform edge length (structured mesh: N = sqrt(Nodes)-1)
    q1 = @(y) -pi*cos(pi/2*y);   % flux on x=+1
    q2 = @(y)  pi*cos(pi/2*y);   % flux on x=-1
    xN = coord(NeuNod,1);
    qfun = @(y) q1(y);           % placeholder, overwritten node-by-node below

    b = zeros(Nodes,1);
    for el = 1:Nelem
        nds = triang(el,:);
        contrib = Area(el)/3;
        for k = 1:3
            xk = coord(nds(k),1); yk = coord(nds(k),2);
            b(nds(k)) = b(nds(k)) + contrib*force(xk,yk);
        end
    end
    % Neumann contribution, split by boundary side (x=+1 vs x=-1)
    isRight = abs(coord(NeuNod,1)-1) < 1e-12;
    isLeft  = abs(coord(NeuNod,1)+1) < 1e-12;
    b(NeuNod(isRight)) = b(NeuNod(isRight)) + q1(coord(NeuNod(isRight),2))*h;
    b(NeuNod(isLeft))  = b(NeuNod(isLeft))  + q2(coord(NeuNod(isLeft),2))*h;

    % ---- impose Dirichlet BC --------------------------------------------
    [sysMat,rhs] = imposeBC(Nodes,NDir,DirNod,DirVal,Stiff,b);

    % ---- solve the linear system (Matlab backslash) ----------------------
    uh = sysMat\rhs;

    % ---- error norm -------------------------------------------------------
    [err(lev), Ai] = computeError(Nelem,Nodes,triang,coord,Area,uh,ureal);
    hval(lev) = h;
    nodes_v(lev) = Nodes;

    fprintf('mesh%d : N=%3d  h=%7.5f  nodes=%6d  eps_h=%12.6e\n', ...
        lev-1, sqrt(Nodes)-1, h, Nodes, err(lev));

    if lev == plot_level
        plot_data.coord  = coord;
        plot_data.triang = triang;
        plot_data.uh     = uh;
        plot_data.uex    = ureal(coord(:,1),coord(:,2));
        plot_data.Nodes  = Nodes;
    end
end

% ---- convergence table --------------------------------------------------
fprintf('\n%-6s%-8s%-10s%-16s%-10s%-10s\n','lev','N','h','eps_h','ratio','p_exp');
ratio = nan(nlev,1); p_exp = nan(nlev,1);
for lev = 2:nlev
    ratio(lev) = err(lev-1)/err(lev);
    p_exp(lev) = log2(ratio(lev));
end
for lev = 1:nlev
    fprintf('%-6d%-8d%-10.5f%-16.6e%-10.3f%-10.3f\n', ...
        lev-1, sqrt(nodes_v(lev))-1, hval(lev), err(lev), ratio(lev), p_exp(lev));
end

% ---- convergence plot (log-log) ------------------------------------------
figure('Name','FEM convergence');
loglog(hval,err,'o-','LineWidth',2,'MarkerSize',7); hold on;
Cref = err(1)/hval(1)^2;
loglog(hval, Cref*hval.^2, 'k--');
xlabel('mesh size h'); ylabel('L^2 error \epsilon_h');
legend('\epsilon_h (FEM)','slope 2 (reference)','Location','northwest');
title('FEM convergence (P1 elements)'); grid on;
print(gcf, 'fig_convergence_loglog.png', '-dpng', '-r150');

% ---- solution / error surface + contour plots (mesh2) ----------------------
figure('Name','Numerical vs exact solution');
subplot(1,3,1);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uh);
shading interp; title('u_h (numerical)'); xlabel('x'); ylabel('y');
subplot(1,3,2);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uex);
shading interp; title('u (exact)'); xlabel('x'); ylabel('y');
subplot(1,3,3);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uh-plot_data.uex);
shading interp; title('error u_h - u'); xlabel('x'); ylabel('y');
print(gcf, 'fig_solution_surfaces.png', '-dpng', '-r150');

figure('Name','Contour plots');
subplot(1,3,1);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uh);
view(2); shading interp; colorbar; title('u_h'); axis equal tight;
subplot(1,3,2);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uex);
view(2); shading interp; colorbar; title('u'); axis equal tight;
subplot(1,3,3);
trisurf(plot_data.triang, plot_data.coord(:,1), plot_data.coord(:,2), plot_data.uh-plot_data.uex);
view(2); shading interp; colorbar; title('u_h - u'); axis equal tight;
print(gcf, 'fig_contours.png', '-dpng', '-r150');

% ---- PCG convergence study (Sec 3, optional point) -------------------------
figure('Name','PCG convergence'); hold on;
cmap = jet(nlev);
for lev = 1:nlev
    meshdir = meshdirs{lev};
    [Nelem,Nodes,triang,coord,NDir,DirNod,DirVal,NNeu,NeuNod,ureal,force] = ...
        input_data(meshdir);
    [Bloc,Cloc,Area] = localBasis(Nelem,triang,coord);
    Stiff = stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area);
    h = 2/(sqrt(Nodes)-1);
    q1 = @(y) -pi*cos(pi/2*y); q2 = @(y) pi*cos(pi/2*y);
    b = zeros(Nodes,1);
    for el = 1:Nelem
        nds = triang(el,:); contrib = Area(el)/3;
        for k = 1:3
            xk = coord(nds(k),1); yk = coord(nds(k),2);
            b(nds(k)) = b(nds(k)) + contrib*force(xk,yk);
        end
    end
    isRight = abs(coord(NeuNod,1)-1) < 1e-12;
    isLeft  = abs(coord(NeuNod,1)+1) < 1e-12;
    b(NeuNod(isRight)) = b(NeuNod(isRight)) + q1(coord(NeuNod(isRight),2))*h;
    b(NeuNod(isLeft))  = b(NeuNod(isLeft))  + q2(coord(NeuNod(isLeft),2))*h;
    [sysMat,rhs] = imposeBC(Nodes,NDir,DirNod,DirVal,Stiff,b);

    L = ichol0(sysMat);
    resvec = [];
    [~,~,~,~,resv] = pcg(sysMat, rhs, 1e-10, 3000, L, L');
    semilogy(1:numel(resv), resv/norm(rhs), 'Color', cmap(lev,:), 'LineWidth',1.5, ...
        'DisplayName', sprintf('mesh%d', lev-1));
end
xlabel('PCG iteration k'); ylabel('||r_k||/||b||');
title('PCG convergence (IC(0) preconditioner)');
legend show; grid on; set(gca,'YScale','log');
print(gcf, 'fig_pcg_convergence.png', '-dpng', '-r150');

disp('Done.');
