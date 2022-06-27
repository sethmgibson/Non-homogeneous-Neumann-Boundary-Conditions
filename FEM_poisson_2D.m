function FEM_poisson_2D()
% Code edited from:
% https://github.com/weizhanghuang/MMPDElab/blob/master/examples/ex2d_heat.m
if (~isdeployed)
  addpath('../src_MMPDElab');
end
% set the basic parameters
   jmax = 31;
   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   nn = 10;
   t = 0;
   tf = 20; % tf = 50; for better convergence
   dt0 = 1e-3;
   dtmax = 0.1;
% set the initial meshes, find the indices of the corner points and fix them
   kmax = jmax;
   [X,tri] = MovMesh_rect2tri(linspace(0,1,jmax), linspace(0,1,kmax), 1);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);
   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;
   % find the indices of the corner points and fix them
   corners = [0, 0; 1, 0; 1, 1; 0, 1];
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
% set initial conditions and compute the initial adjusted mesh
   % set the initial solution
   U = zeros(Nv,npde);   
   % generate initial adjusted mesh
   if (moving_mesh)   
      for n=1:nn  
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([0,1],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
         X = Xnew;     
         figure(1)
         triplot(tri,X(:,1),X(:,2),'Color','r')
         axis([0 1 0 1]);  
         axis square;
         drawnow;
      end
   end  
% define PDE system and BCs   
   % define boundary types
   pdedef.bftype = zeros(Nbf,npde);           
   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
% perform integration (MP)      
   dt = dt0;
   DT = zeros(20000,2);
   n = 0;   
   tcpu = cputime;   
   while true
      % move the mesh      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([t,t+dt],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      else
         Xnew = X;
      end
      Xdot = (Xnew-X)/dt;      
      % integrate physical PDEs      
      [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef);      
      % update
      X = X + dt0*Xdot;
      U = Unew;
      n = n + 1;     
      DT(n,:) = [t, dt0];      
      t = t + dt0;
      dt = min(dtmax,dt1);
      if (t+dt>tf), dt=tf-t; end      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e', ...
                     n,t,dt0,dt1);      
      figure(2)
      clf
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([0 1 0 1]);
      axis square;
      drawnow;     
      if (t>=tf-100*eps || dt < 100*eps), break; end      
   end   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);
% output
   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1))   
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2));   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);       
end

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)
    F = 10.*exp( -50*( (x(:,1)-0.5).^2 + (x(:,2)-0.5).^2 ) );
    k = 1; % k = 100; for better convergence
    F = -k*ut(:,1).*v(:) - du(:,1).*dv(:,1) - du(:,2).*dv(:,2) - F.*v(:); 
end

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)
   G = -sin(5*x(:,1)).*v(:);
end