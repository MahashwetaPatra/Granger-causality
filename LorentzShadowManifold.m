%
% NOTE: This is for calculating orbit for t=100 for one sigle realization in Lorentz system s7
%       Shadow manifold M_x, M_y is plotted with the coordinate (x(t), x(t-\tau) and x(t-2\tau)) and (y(t), y(t-\tau) and y(t-2\tau)) 
% HIST:  - June, 2021: Created by Patra
%=========================================================================
clc; close all; clear all;
%initializing x,y,z
x(1)=-0.01; y(1)=0.01; z(1)=0.01;
dt = 0.001;
h=0.001;   %step size
t=0:h:100;
t1=0.1:h:100.1; % new time series with time lag 0.1
t2=0.2:h:100.2; % new time series with time lag 0.2
f1=@(t,x,y,z) 10*(-x+y);  %ode 
g1=@(t,x,y,z) 28*x-y-x*z; 
p1=@(t,x,y,z) x*y-(8/3)*z;
n=(length(t)-1);
for i=1:n %loop
    y0=y(i);
    k1=f1(t(i),x(i),y(i),z(i));
    l1=g1(t(i),x(i),y(i),z(i));
    m1=p1(t(i),x(i),y(i),z(i));
    k2=f1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
    l2=g1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
    m2=p1(t(i)+h/2,(x(i)+0.5*k1*h),(y(i)+(0.5*l1*h)),(z(i)+(0.5*m1*h)));
    k3=f1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
    l3=g1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
    m3=p1(t(i)+h/2,(x(i)+0.5*k2*h),(y(i)+(0.5*l2*h)),(z(i)+(0.5*m2*h)));
    k4=f1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
    l4=g1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
    m4=p1(t(i)+h,(x(i)+k3*h),(y(i)+l3*h),(z(i)+m3*h));
    x(i+1) = x(i) + h*(k1 +2*k2  +2*k3   +k4)/6; %final equations
    y(i+1) = y(i) + h*(l1  +2*l2   +2*l3    +l4)/6;
    z(i+1) = z(i) + h*(m1+2*m2 +2*m3  +m4)/6;        
end
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontSize',14,'LineWidth',2.75)
figure(1);
plot3(x,y,z, '.','markersize',5)
hold on;
plot3(x(9000), y(9000), z(9000), '.b', 'markersize', 20)
xlabel('normalized:n tangential wind u');
ylabel('normalized vertical wind v');
zlabel('normalized temperature anomaly b')
grid on
view(58,9)
figure(2);
k=201:1:n;
view(58,9)
plot3(x(k),x(k-100),x(k-200), '.','markersize',5)
hold on;
plot3(x(9000), x(8900), x(8800), '.r', 'markersize', 20)
figure(3);
k=201:1:n;
view(58,9)
plot3(y(k),y(k-100),y(k-200), '.','markersize',5)
hold on;
plot3(y(9000), y(8900), y(8800), '.g', 'markersize', 20)