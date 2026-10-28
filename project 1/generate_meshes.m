function generate_meshes(outdir)
% GENERATE_MESHES  Build a nested sequence of 5 uniformly refined
% triangular meshes of the square Omega=[-1,1]x[-1,1], following the
% refinement scheme of the project (each refinement level
% connects edge midpoints; for a structured criss-cross grid this is
% equivalent to doubling the number of subdivisions per side N).
%
% Writes, for each level k=0..4, a folder meshK/ containing:
%   triang.dat  - element/node connectivity (Nel x 4, last col = zone id)
%   xy.dat      - node coordinates (Nodes x 2)
%   dirtot.dat  - indices of ALL boundary nodes
%   dirnod.dat  - indices of Dirichlet boundary nodes (y = +-1)
%   neunod.dat  - indices of Neumann boundary nodes   (x = +-1, y~=+-1)
%
% Usage: generate_meshes('.')

if nargin < 1, outdir = '.'; end
Nlist = [4 8 16 32 64];

for lev = 0:4
    N = Nlist(lev+1);
    x = linspace(-1,1,N+1);
    y = linspace(-1,1,N+1);
    [X,Y] = meshgrid(x,y);
    coord = [X(:), Y(:)];

    nidx = @(i,j) j*(N+1) + i + 1;  % 1-based linear index, i,j in [0,N]

    triang = zeros(2*N*N,4);
    e = 0;
    for j = 0:N-1
        for i = 0:N-1
            n00 = nidx(i,   j);
            n10 = nidx(i+1, j);
            n11 = nidx(i+1, j+1);
            n01 = nidx(i,   j+1);
            e = e+1; triang(e,:) = [n00 n10 n11 1];
            e = e+1; triang(e,:) = [n00 n11 n01 1];
        end
    end

    tol = 1e-12;
    xs = coord(:,1); ys = coord(:,2);
    dirnod = find(abs(ys-1)<tol | abs(ys+1)<tol);
    neunod = find( (abs(xs-1)<tol | abs(xs+1)<tol) & abs(ys-1)>tol & abs(ys+1)>tol );
    dirtot = sort([dirnod; neunod]);

    meshdir = fullfile(outdir, sprintf('mesh%d',lev));
    if ~exist(meshdir,'dir'), mkdir(meshdir); end
    write_dat(fullfile(meshdir,'triang.dat'), triang);
    write_dat(fullfile(meshdir,'xy.dat'),     coord);
    write_dat(fullfile(meshdir,'dirtot.dat'), dirtot);
    write_dat(fullfile(meshdir,'dirnod.dat'), dirnod);
    write_dat(fullfile(meshdir,'neunod.dat'), neunod);

    fprintf('mesh%d : N=%2d  nodes=%5d  elements=%5d\n', lev, N, size(coord,1), size(triang,1));
end
end

function write_dat(filename, M)
% Portable ASCII matrix writer (works identically in Matlab and Octave,
% unlike writematrix which is Matlab-only / recent versions only).
fid = fopen(filename,'w');
fmt = [repmat('%.10g ', 1, size(M,2)) '\n'];
for r = 1:size(M,1)
    fprintf(fid, fmt, M(r,:));
end
fclose(fid);
end
