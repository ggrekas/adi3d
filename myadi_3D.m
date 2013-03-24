% Alternative Direction Implicit (ADI) method implementation for solving
% parabolic partial differential equations, with Neumann boundary
% conditions:
% u_t = a*(u_xx + u_yy + u_zz) + Cg*\Nabla[u \Nabla(g)]+Cphi*\Nabla[u \Nabla(phi)]
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
function myadi_3D(u, a, C, f, h_t, Cg, g, Cphi, phi)
N = size(u,1);
h = 1/(N-1);

f = h*h*h_t*f;
% applies operators Ax, Ay, Az respectively
[sub_diag_x, diag_x, hyp_diag_x] = compute_x_diags(a, Cg, g, Cphi, phi, C, h); %TODO remove 1/h^2 in derivatives
[sub_diag_y, diag_y, hyp_diag_y] = compute_y_diags(a, Cg, g, Cphi, phi, C, h);
[sub_diag_z, diag_z, hyp_diag_z] = compute_z_diags(a, Cg, g, Cphi, phi, C, h);

sub_diag_x = 0.5*h_t*sub_diag_x;
hyp_diag_x = 0.5*h_t*hyp_diag_x;
diag_x = 0.5*h_t*diag_x;

sub_diag_y = h_t*sub_diag_y;
hyp_diag_y = h_t*hyp_diag_y;
diag_y = h_t*diag_y;

sub_diag_z = h_t*sub_diag_y;
hyp_diag_z = h_t*hyp_diag_y;
diag_z = h_t*diag_y;


u = x_sweep(u, f, sub_diag_x, diag_x, hyp_diag_x, sub_diag_y, diag_y, hyp_diag_y, ...
   sub_diag_z, diag_z, hyp_diag_z, h_t);

return;

function u = x_sweep(u, f, sub_diag_x, diag_x, hyp_diag_x, sub_diag_y, diag_y, hyp_diag_y, ...
   sub_diag_z, diag_z, hyp_diag_z, h_t)

rhs = rhs_calculation_x(u, f, sub_diag_x, diag_x, hyp_diag_x, sub_diag_y, diag_y, hyp_diag_y, ...
   sub_diag_z, diag_z, hyp_diag_z, h_t);




return;

function rhs = rhs_calculation_x(u, f, sub_diag_x, diag_x, hyp_diag_x, sub_diag_y, diag_y, hyp_diag_y, ...
   sub_diag_z, diag_z, hyp_diag_z, h_t)
N = size(u,1);
h= 1/(N-1);
rhs = zeros( size(u) );

for i=1:N
   for j=1:N
      for k=1:N
         rhs(i, j, k) = 
      end
   end
end

rhs(2:end-1, 2:end-1, 2:end-1) = (h*h + sub_diag_x(1:end-2, 1:end-2, 1:end-2) +...
   hyp_diag)

return;



function [sub_diag, diag, hyp_diag] = compute_x_diags(a, Cg, g, Cphi, phi,...
   C, h)
a_x = derivative_x(a);
g_x = derivative_x(g);
phi_x = derivative_x(phi);

u_x_coeff = 0.5*h*(a_x + Cg.*g_x + Cphi.*phi_x); 
hyp_diag = a + u_x_coeff;
sub_diag = a - u_x_coeff;

diag = -2*a + h*h*( Cg.*derivative_xx(g) + Cphi.*derivative_xx(phi) +...
   1/3*C);

return;

function [sub_diag, diag, hyp_diag] = compute_y_diags(a, Cg, g, Cphi, phi,...
   C, h)
a_y = derivative_y(a);
g_y = derivative_y(g);
phi_y = derivative_y(phi);

u_y_coeff = 0.5*h*(a_y + Cg.*g_y + Cphi.*phi_y); 
hyp_diag = a + u_y_coeff;
sub_diag = a - u_y_coeff;

diag = -2*a + h*h*( Cg.*derivative_yy(g) + Cphi.*derivative_yy(phi) +...
   1/3*C);

return;

function [sub_diag, diag, hyp_diag] = compute_z_diags(a, Cg, g, Cphi, phi,...
   C, h)
a_z = derivative_z(a);
g_z = derivative_z(g);
phi_z = derivative_z(phi);

u_z_coeff = 0.5*h*(a_z + Cg.*g_z + Cphi.*phi_z); 
hyp_diag = a + u_z_coeff;
sub_diag = a - u_z_coeff;

diag = -2*a + h*h*( Cg.*derivative_zz(g) + Cphi.*derivative_zz(phi) +...
   1/3*C);

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
