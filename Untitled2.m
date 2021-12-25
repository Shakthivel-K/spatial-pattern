% % defining functions:
% global
a=12.0;
b=0.2;
d=1.0;
sigma=50.0;
lp=0.47;% 0.43
R=4.0;
lf=R.*lp;
W=3;
h=0.05;
syms x y;
syms u v;
syms U X;
w=@(x,y) (W.*((sinpi(2.*x./lf)+sinpi((2./lf).*(0.5.*x+3.^0.5.*y./2)+1./6)+sinpi((2./lf).*(0.5.*x-3.^0.5.*y./2)+1./3)).*(2./9)+(1./3)));
f=@(t,u,v,x,y) (a-u-(4.*u.*v./(1+u.^2))-w(x,y));
g=@(t,u,v,x,y) (sigma.*b.*(u-(u.*v./(1+u.^2))+w(x,y)));
F=@(t,U,X) ([f(t,U(1),U(2),X(1),X(2))+laplacian(U(1),X);g(t,U(1),U(2),X(1),X(2))+d.*laplacian(U(2),X)]);
%critical point 
u0=a./5-W;
v0=a./5+(1+u0.^2)./u0;
zeta=[u-u0;v-v0];
syms Z;
[T,Z] = ode45(F, [0 100], zeta)



%boundary condition 
%x=0

