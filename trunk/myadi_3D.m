% Alternative Direction Implicit (ADI) method implementation for solving
% parabolic partial differential equations, with Neumann boundary
% conditions:
% u_t = a*(u_xx + u_yy + u_zz) + a2*\Nabla[u \Nabla(g2)]+a3*\Nabla[u \Nabla(g3)]
% + c*u + f
%
% function myadi(u, a, c, f, h_t, a2, g2, a3, g3)
%
% 
%
% Input Arguments:
%	u:    [NxNxN]: the input function. u(x,y,z, n)
%	a:  	[1x1] or [NxNxN]: the diffusion coefficient.
%	a2:  	[1x1] or [NxNxN].
%	g2:  	[NxNxN].
%	a3:  	[1x1] or [NxNxN].
%	g3:  	[NxNxN].
%	C:  	[NxNxN]. c(x,y,z,n+1/2)
%	f: [NxN]: f(x,y,n+1/2):=0.5*[ f(x,y,n)+f(x,y,n+1) ]
%
% Output Arguments:
%	u:	[NxNxN]: the function u(x,y,z,n+1)
%
function myadi_3D
x= zeros(5,5,2);

for j=1:5
   for k=1:2
      x(1:end,j,k)=j;
   end
end

dx= derivative_x(x, 'no');
dy= derivative_y(x, 'no');
return;


function df = derivative_xx(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(2:end-1,:,:)= (f(3:end,:,:) -2*df(2:end-1,:,:)+ f(1:end-2,:,:))*(N-1)^2; 
df(1,:,:)= 2*( f(2,:,:) -f(1,:,:) )* (N-1)^2;
df(end,:,:)= 2*( f(end-1,:,:) -f(end,:,:) )* (N-1)^2;


% if( strcmp(boundaryCond,'dirichlet') )
%    df(1,:,:)= ( f(2,:,:) - f(1,:,:) )* (N-1);
%    df(end,:,:)= ( f(end,:,:) - f(end-1,:,:) )* (N-1);
% end
return;


function df = derivative_yy(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1,:)= (f(:,3:end,:) -2*df(:,2:end-1,:)+ f(:,1:end-2,:))*(N-1)^2; 
df(:,1,:)= 2*( f(:,2,:) -f(:,1,:) )* (N-1)^2;
df(:,end,:)= 2*( f(:,end-1,:) -f(:,end,:) )* (N-1)^2;

return;


function df = derivative_zz(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,:,2:end-1)= (f(:,:,3:end) -2*df(:,:,2:end-1)+ f(:,:,1:end-2))*(N-1)^2; 
df(:,:,1)= 2*( f(:,:,2) -f(:,:,1) )* (N-1)^2;
df(:,:,end)= 2*( f(:,:,end-1) -f(:,:,end) )* (N-1)^2;

return;

function df = derivative_x(f, boundaryCond)
N= size(f,1); %h = 1/(N-1)


df= zeros(size(f));

df(2:end-1,:,:)= (f(3:end,:,:) - f(1:end-2,:,:))*0.5*(N-1); 
if( strcmp(boundaryCond,'dirichlet') )
   df(1,:,:)= ( f(2,:,:) - f(1,:,:) )* (N-1);
   df(end,:,:)= ( f(end,:,:) - f(end-1,:,:) )* (N-1);
end
return;

function df = derivative_y(f, boundaryCond)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1,:)= (f(:,3:end,:) - f(:,1:end-2,:))*0.5*(N-1); 
if( strcmp(boundaryCond,'dirichlet') )
   df(:,1,:)= ( f(:,2,:) - f(:,1,:) )* (N-1);
   df(:,end,:)= ( f(:,end,:) - f(:,end-1,:) )* (N-1);
end
%df= derivative_x(permute(f, [2,1,3]), boundaryCond);
return;

function df = derivative_z(f, boundaryCond)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,:,2:end-1)= (f(:,:,3:end) - f(:,:,1:end-2))*0.5*(N-1); 
if( strcmp(boundaryCond,'dirichlet') )
   df(:,:,1)= ( f(:,:,2) - f(:,:,1) )* (N-1);
   df(:,:,end)= ( f(:,:,end) - f(:,:,end-1) )* (N-1);
end
%df= derivative_x(permute(f, [3, 1, 2]), boundaryCond);
return;
