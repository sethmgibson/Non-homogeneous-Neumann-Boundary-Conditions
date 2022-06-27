% MAT 494 - Introduction to Numerical Methods for Partial Differential Equations
% Based on Example 5.4.1: Finite volume method by Xianping Li
% Using FVTool to solve a 2D steady-state diffusion equation on \Omega=(0,1)x(0,1)
% https://github.com/simulkade/FVTool

clc; clear;

% add paths
addpath('FVTool');
FVToolStartUp

% Define the domain and create a mesh structure
L = 1;  % domain length
Nx = 2^6; % number of cells
dx = L/Nx;
m = createMesh2D(Nx,Nx, L,L);

x = m.cellcenters.x;
y = m.cellcenters.y;

% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = -1; BC.left.b(:)=0; BC.left.c(:)=bc(0,y); % left boundary
BC.right.a(:) = 1; BC.right.b(:)=0; BC.right.c(:)= bc(1,y); % right boundary
BC.top.a(:) = 1; BC.top.b(:)=0; BC.top.c(:)=bc(x,1); % top boundary
BC.bottom.a(:) = -1; BC.bottom.b(:)=0; BC.bottom.c(:)=bc(x,0); % bottom boundary

% define the transfer coeffs
D_x = -1;
D_y = -1;
D = createFaceVariable(m, [D_x, D_y]);
Mdiff = diffusionTerm(D);

% define source term
[X, Y] = ndgrid(x, y);
S = @(X,Y)source(X,Y);
s1 = constantSourceTerm(createCellVariable(m,S(X,Y)));

[Mbc, RHSbc] = boundaryCondition(BC);

% Subtracting mean (works better for higher number of cells).
M = (Mdiff+Mbc) - mean(Mdiff+Mbc); 
RHS = (-s1+RHSbc) - mean(-s1+RHSbc);
c = solvePDEpi(m,M, RHS);

figure(1);visualizeCells(c);caxis([0,10]);drawnow;


% view 3D plot
% add ghost points to x and y
x = [-0.5*dx; x; L+0.5*dx];
y = [-0.5*dx; y; L+0.5*dx];
[X, Y] = ndgrid(x, y);

% fix corner points as ghost point
u = c.value;
u(1,1) = u(1,2);
u(1,end) = u(1,end-1);
u(end,1) = u(end,2);
u(end,end) = u(end,end-1);

figure(2); surf(X, Y, u);
%axis([0 L 0 L 0 10]);


function f = source(x,y)
    f = 10*exp(-50*((x-0.5).^2+(y-0.5).^2));
end

function g = bc(x,y)
    g = -sin(5*x);
end

function phi = solvePDEpi(MS, M, RHS, varargin)
%SOLVEPDE solves the linear system M x \phi = RHS and returns the value of
%\phi reshaped based on the structure of the MeshStructure variable. The
%default solver is the matlab '\' linear solver

% Written by Ali A. Eftekhari
% See the license file

n = MS.dimension;
N = MS.dims;

i = floor(0.5*N);
j = floor(0.5*N);
        
x(j, i) = 0.0;
x = M\RHS;

% improved fix point
%         P = M/(M'*M)\M';
%         f1 = P*RHS;
%         f1-RHS;
%         x = M\f1;

if (n>=2)
    phival = reshape(x, N+2);
else
    phival = reshape(x, [N(1)+2 1]);
end

phi=CellVariable(MS, phival);

end