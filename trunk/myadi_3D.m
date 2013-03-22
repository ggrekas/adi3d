% Alternative Direction Implicit (ADI) method implementation for solving
% parabolic partial differential equations, with Neumann boundary
% conditions:
% u_t = a*(u_xx + u_yy + u_zz) + Cg*\Nabla[u \Nabla(g)]+Ch*\Nabla[u \Nabla(h)]
% + c*u + f
%
% function myadi(u, a, c, f, h_t, Cg, g, Ch, h)
%
% 
%
% Input Arguments:
%	u:    [NxNxN]: the input function. u(x,y,z, n)
%	a:  	[1x1] or [NxNxN]: the diffusion coefficient.
%	Cg:  	[1x1] or [NxNxN].
%	g:  	[NxNxN].
%	Ch:  	[1x1] or [NxNxN].
%	h:  	[NxNxN].
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

dx= derivative_x(x);
dy= derivative_y(x);

[sub_diag_x, diag_x, hyp_diag_x] = compute_hyp_sub_diag_x(a, Cg, g, Ch, h, h_t);
return;


function [sub_diag, diag, hyp_diag] = compute_hyp_diag_x(a, Cg, g, Ch, h, h_t)
g_x = derivative_x(g);
h_x = derivative_x(h);

hyp_diag = 0.5*ht*()



return;



function df = derivative_xx(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(2:end-1,:,:)= ( f(3:end,:,:) -2*f(2:end-1,:,:)+ f(1:end-2,:,:) ); 
df(1,:,:)= 2*( f(2,:,:) -f(1,:,:) );
df(end,:,:)= 2*( f(end-1,:,:) -f(end,:,:) );

df = (N-1)^2*df;


% if( strcmp(boundaryCond,'dirichlet') )
%    df(1,:,:)= ( f(2,:,:) - f(1,:,:) )* (N-1);
%    df(end,:,:)= ( f(end,:,:) - f(end-1,:,:) )* (N-1);
% end
return;


function df = derivative_yy(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1,:)= f(:,3:end,:) -2*f(:,2:end-1,:)+ f(:,1:end-2,:); 
df(:,1,:)= 2*( f(:,2,:) -f(:,1,:) );
df(:,end,:)= 2*( f(:,end-1,:) -f(:,end,:) ) ;

df = (N-1)^2*df;
return;


function df = derivative_zz(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,:,2:end-1)= f(:,:,3:end) -2*df(:,:,2:end-1)+ f(:,:,1:end-2); 
df(:,:,1)= 2*( f(:,:,2) -f(:,:,1) );
df(:,:,end)= 2*( f(:,:,end-1) -f(:,:,end) );

df = (N-1)^2*df;
return;

function df = derivative_x(f)
N= size(f,1); %h = 1/(N-1)


df= zeros(size(f));

df(2:end-1,:,:)= (f(3:end,:,:) - f(1:end-2,:,:))*0.5*(N-1); 
return;

function df = derivative_y(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,2:end-1,:)= (f(:,3:end,:) - f(:,1:end-2,:))*0.5*(N-1); 
return;

function df = derivative_z(f)
N= size(f,1); %h = 1/(N-1)
df= zeros(size(f));

df(:,:,2:end-1)= (f(:,:,3:end) - f(:,:,1:end-2))*0.5*(N-1); 
return;
