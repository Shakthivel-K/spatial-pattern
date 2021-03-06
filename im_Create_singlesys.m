% clear all MATLAB variables, command windows, and figure
clear all; clc; close all;

% spatial discretization
xright=10; nx=50;
yright=10; ny=floor(nx*yright/xright); 
h=xright/nx;
x=linspace(-0.5*h,xright+0.5*h,nx+2)';
y=linspace(-0.5*h,yright+0.5*h,ny+2);

% set the parameters which govern equation
Du=1;       %diffusion coefficients
Dv=1;    %diffusion coefficients
k1=0.2;      %same as parameters a
k2=12;      %same as parrameters b
ubar=1+0.04*k2^2;
vbar=0.2*k2;
w=1.9;
sigma=50;
% set the parameters (time discretization)
dt=0.1*h^2;
maxit=80000;
nn=maxit;

for pit=1:9    % the number of images
% set the initial condition
rng(pit);
u=ubar+0.1*(2*rand(nx+2,ny+2)-1);
v=vbar+0.1*(2*rand(nx+2,ny+2)-1);
nu=u; nv=v;

% nx=5;
% u1=u0+0.1*(2*rand(nx+2)-1);
% v1=v0+0.1*(2*rand(nx+2)-1);
% zhiu=u1; zhiv=v1;

% numerical scheme
    for it=1:maxit
    % periodic boundary condition
%         u(2:end-1,1)=u(2:end-1,end-1);
%         u(2:end-1,end)=u(2:end-1,2);
%         u(1,:)=u(end-1,:);
%         u(end,:)=u(2,:);
%         v(2:end-1,1)=v(2:end-1,end-1);
%         v(2:end-1,end)=v(2:end-1,2);
%         v(1,:)=v(end-1,:);
%         v(end,:)=v(2,:);
    
    % utilda=@(t,x) ustar*exp(sig*t)*sin(alpha*x);% Substitute in dzhi_udt=d1*lap(zhi_u)+fu*zhi_u+fv*zhi_v;      at (u0,v0)
    % vtilda=@(t,x) vstar*exp(sig*t)*sin(alpha*x);%               dzhi_vdt=sigma*d2*lap(zhi_v)+gu*zhi_u+gv*zhi_v;
    % diffu=(d1*laplacian(u1,[t x]))%(u1,h));
    % diffv=sigma*(d2*laplacian(v1,[t x]))%,h)));
    
    % set the source terms
        F=u(2:end-1,2:end-1).*v(2:end-1,2:end-1) ...
            ./(1+v(2:end-1,2:end-1).^2);
        f=sigma*k1*(v(2:end-1,2:end-1)-F);
        g=k2-v(2:end-1,2:end-1)-4*F;

    % solve the equations
        nu(2:end-1,2:end-1)=u(2:end-1,2:end-1)+dt*(f+Du*sigma*lap(u,h)+w);
        nv(2:end-1,2:end-1)=v(2:end-1,2:end-1)+dt*(g+Dv*lap(v,h)-w);
    
    % solve the equations
    % zhiu(2:end-1,2:end-1)=u1(2:end-1,2:end-1)+dt*(diffu+(A0(1,1).*u1(2:end-1,2:end-1))+(A0(1,2).*v1(2:end-1,2:end-1)));
    % zhiv(2:end-1,2:end-1)=v1(2:end-1,2:end-1)+dt*(diffv+(A0(2,1).*u1(2:end-1,2:end-1))+(A0(2,2).*v1(2:end-1,2:end-1)));
    
    % reset the variables for next step
        u=nu;
        v=nv;
    
    % visualization
        figure(pit);
        
        %print('-djpeg',sprintf('0/pattern_%d',pit));    % storage path
    end
    % % display content
    surf(x(2:end-1),y(2:end-1),u(2:end-1,2:end-1)','linestyle','none');
    axis image;
    view(2);
    set(gca, 'xtick',[], 'ytick',[]);
    box on;
    shading interp;
    drawnow;
    
end
