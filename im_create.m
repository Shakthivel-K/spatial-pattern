% clear all MATLAB variables, command windows, and figure
clear all; clc; close all;

% spatial discretization
xright=10; nx=50;
yright=10; ny=floor(nx*yright/xright); 
h=xright/nx;
x=linspace(-0.5*h,xright+0.5*h,nx+2)';
y=linspace(-0.5*h,yright+0.5*h,ny+2);
alphaa=0.1;
betaa=0.1;
% alphaa=0.5; changing to these values coverts to labrinth
% betaa=0.5;
w=0;
% set the parameters which govern equation
Du2=1.404;       %diffusion coefficients
Dv2=2.088;
Du1=10.719;
Dv1=22.796;
%diffusion coefficients
k1=1.5;      %same as parameters a
k2=18;      %same as parrameters b
u1bar=1+0.04*k2^2;
v1bar=0.2*k2;
u2bar=1+0.04*k2^2;
v2bar=0.2*k2;
sigma=9.0;
% set the parameters (time discretization)
dt=0.1*h^2;
maxit=2;
nn=maxit;

for pit=1:8    % the number of images
% set the initial condition
rng(pit);
u1=u1bar+0.1*(2*rand(nx+2,ny+2)-1);
v1=v1bar+0.1*(2*rand(nx+2,ny+2)-1);
u2=u2bar+0.1*(2*rand(nx+2,ny+2)-1);
v2=v2bar+0.1*(2*rand(nx+2,ny+2)-1);
nu1=u1; nv1=v1;
nu2=u2; nv2=v2;

% nx=5;
% u1=u0+0.1*(2*rand(nx+2)-1);
% v1=v0+0.1*(2*rand(nx+2)-1);
% zhiu=u1; zhiv=v1;

% numerical scheme
    for it=1:maxit
    % periodic boundary condition
%         u1(2:end-1,1)=u1(2:end-1,end-1);
%         u1(2:end-1,end)=u1(2:end-1,2);
%         u1(1,:)=u1(end-1,:);
%         u1(end,:)=u1(2,:);
%         v1(2:end-1,1)=v1(2:end-1,end-1);
%         v1(2:end-1,end)=v1(2:end-1,2);
%         v1(1,:)=v1(end-1,:);
%         v1(end,:)=v1(2,:);
%         u2(2:end-1,1)=u2(2:end-1,end-1);
%         u2(2:end-1,end)=u2(2:end-1,2);
%         u2(1,:)=u2(end-1,:);
%         u2(end,:)=u2(2,:);
%         v2(2:end-1,1)=v2(2:end-1,end-1);
%         v2(2:end-1,end)=v2(2:end-1,2);
%         v2(1,:)=v2(end-1,:);
%         v2(end,:)=v2(2,:);
    
    % utilda=@(t,x) ustar*exp(sig*t)*sin(alpha*x);% Substitute in dzhi_udt=d1*lap(zhi_u)+fu*zhi_u+fv*zhi_v;      at (u0,v0)
    % vtilda=@(t,x) vstar*exp(sig*t)*sin(alpha*x);%               dzhi_vdt=sigma*d2*lap(zhi_v)+gu*zhi_u+gv*zhi_v;
    % diffu=(d1*laplacian(u1,[t x]))%(u1,h));
    % diffv=sigma*(d2*laplacian(v1,[t x]))%,h)));
    
    % set the source terms
        F1=u1(2:end-1,2:end-1).*v1(2:end-1,2:end-1)...
        ./(1+v1(2:end-1,2:end-1).^2);
        F2=u2(2:end-1,2:end-1).*v2(2:end-1,2:end-1)...
        ./(1+v2(2:end-1,2:end-1).^2);
        f1=sigma*k1*(v1(2:end-1,2:end-1)-F1+w)+alphaa.*(u2(2:end-1,2:end-1)-u1(2:end-1,2:end-1));
        g1=k2-v1(2:end-1,2:end-1)-4*F1-w+betaa.*(v2(2:end-1,2:end-1)-v1(2:end-1,2:end-1));
        f2=sigma*k1*(v2(2:end-1,2:end-1)-F2+w)+alphaa.*(u2(2:end-1,2:end-1)-u1(2:end-1,2:end-1));
        g2=k2-v2(2:end-1,2:end-1)-4*F2-w+betaa.*(v2(2:end-1,2:end-1)-v1(2:end-1,2:end-1));
    % solve the equations
        nu1(2:end-1,2:end-1)=u1(2:end-1,2:end-1)+dt*(f1+sigma*Du1*lap(u1,h));
        nv1(2:end-1,2:end-1)=v1(2:end-1,2:end-1)+dt*(g1+Dv1*lap(v1,h));
        nu2(2:end-1,2:end-1)=u2(2:end-1,2:end-1)+dt*(f2+sigma*Du2*lap(u2,h));
        nv2(2:end-1,2:end-1)=v2(2:end-1,2:end-1)+dt*(g2+Dv2*lap(v2,h));
    
    % solve the equations
    % zhiu(2:end-1,2:end-1)=u1(2:end-1,2:end-1)+dt*(diffu+(A0(1,1).*u1(2:end-1,2:end-1))+(A0(1,2).*v1(2:end-1,2:end-1)));
    % zhiv(2:end-1,2:end-1)=v1(2:end-1,2:end-1)+dt*(diffv+(A0(2,1).*u1(2:end-1,2:end-1))+(A0(2,2).*v1(2:end-1,2:end-1)));
    
    % reset the variables for next step
        u1=nu1;
        v1=nv1;
        u2=nu2;
        v2=nv2;
    
    % visualization
        figure(pit);
        
        %print('-djpeg',sprintf('0/pattern_%d',pit));    % storage path
    end
    % % display content
    surf(x(2:end-1),y(2:end-1),u1(2:end-1,2:end-1)','linestyle','none');
    axis image;
    view(2);
    set(gca, 'xtick',[], 'ytick',[]);
    box on;
    shading interp;
    drawnow;
    
end



