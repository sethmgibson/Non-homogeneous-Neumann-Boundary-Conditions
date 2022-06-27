% MAT 494 - Introduction to Numerical Methods for Partial Differential Equations
% Example 4.3.2: Finite difference method
% Xianping Li

M=2^5
N=M

% Solve the 2D poisson equation on space domain (0,1)x(0,1)
% u_{xx} + u_{yy} = f(x,y) = 10*exp(-50*((x-0.5)^2+(y-0.5)^2))
% with Neumann boundary condition du/dn = g = -sin(5x)

% different ways to deal with singular matrix

dx = 1.0/M;
dy = 1.0/N;

a = 1/dx/dx;
b = 1/dy/dy;
c = -2*(a+b);
d = b;
e = a;

% create the nodes in 2D space
x = linspace(0,1,M+1);
x = x'; % convert to a column vector
y = linspace(0,1,N+1);
[X, Y] = meshgrid(x,y);
X;
Y;

% plot the grids with u=0
u = 0*X;
figure(1)
clf;
surf(X,Y,u);
xlabel('x');
ylabel('y');

% initialize the vector u and right-side vector
% use matrix form first then reshape as a vector
f = source(X, Y); 
    % boundary points at x=0;
    f(:,1) = bc(X(:,1),Y(:,1))*dx;
    % boundary points at x=1;
    f(:,M+1) = bc(X(:,M+1),Y(:,M+1))*dx;
    % boundary points at y=0;
    f(1,:) = bc(X(1,:),Y(1,:))*dy;
    % boundary points at y=1;
    f(N+1,:) = bc(X(N+1,:),Y(N+1,:))*dy;

u = reshape(u, [(N+1)*(M+1),1]);

%% initialize the matrix A, including boundary points
A = zeros((N+1)*(M+1), (N+1)*(M+1));
% first block for boundary x=0
i = 1;
for j = 1 : N+1
    k = (N+1)*(i-1) + j;
    A(k, k) = 1.0; 
    A(k, k+(N+1)) = -1.0;
end
% middle blocks
for i = 2 : M
    % bottom boundary
    j = 1;
    k = (N+1)*(i-1) + j;
    A(k, k) = 1.0;
    A(k, k+1) = -1.0;
    
    % top boundary
    j = N+1;
    k = (N+1)*(i-1) + j;
    A(k, k) = 1.0;
    A(k, k-1) = -1.0;
    
    % middle points
    for j = 2 : N
        % diagonal entries
        k = (N+1)*(i-1) + j;
        A(k, k-(N+1)) = a;
        A(k, k-1) = b;
        A(k, k) = c;
        A(k, k+1) = d;
        A(k, k+(N+1)) = e;
    end
end

% last block for boundary x=1
i = M + 1;
for j = 1 : N+1
    k = (N+1)*(i-1) + j;
    A(k, k) = 1.0; 
    A(k, k-(N+1)) = -1.0;
end

%f = f - mean(f);

        % problem: the matrix A is singular due to the pure Neumann boundary condition
        % Method II: elimilating a point u_{i,j} from the governing equation and
        % replace it with a constant value, see u_{i,j}=1
        i = floor(0.5*M);
        j = floor(0.5*N);
        k = (N+1)*(i-1) + j;
        A(k, k-(N+1)) = 0;
        A(k, k-1) = 0;
        A(k, k) = 1.0;
        A(k, k+1) = 0;
        A(k, k+(N+1)) = 0;
        
        cond(A)
        
        f(j, i) = 1.0;
        f = reshape(f, [(N+1)*(M+1),1]);
        
        % improved fix point
        P = A*inv(A'*A)*A';
        f = P*f;
        u = A\f

% reshape u to matrix form
u = reshape(u, [N+1, M+1]);

% show results
surf(X,Y,u);
xlabel('x');
ylabel('y');


function f = source(x,y)

    f = 10*exp(-50*((x-0.5).^2+(y-0.5).^2));

end


function g = bc(x,y)

    g = -sin(5*x);

end