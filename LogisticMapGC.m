%#################### Logistic Map GC ########################
%
% NOTE:     For 2D delayed Logistic map, it checks the time series in different time zones
%           and calculates the correlation coefficient between the variables in those time zones
%          
%         
%=========================================================================

clc;close all;clear all;
x(1) =0.4;y(1) =0.2;
n=1:1:24001;
for j = 1:24000
    x(j+1) = x(j)*(3.8-3.8*x(j)-0.02*y(j));
    y(j+1) = y(j)*(3.5-3.5*y(j)-0.1*x(j));
end
size(n);
size(x);

%plot the time series for different time zones
figure(1);
subplot(1,3,1)
plot(n(377:388),x(377:388),'-b','linewidth',2)
hold on;
plot(n(377:388),y(377:388),'-r','linewidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
grid on
subplot(1,3,2)
plot(n(903:920),x(903:920),'-b','linewidth',2)
hold on;
plot(n(903:920),y(903:920),'-r','linewidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
grid on
subplot(1,3,3)
plot(n(415:440),x(415:440),'-b','linewidth',2)
hold on;
plot(n(415:440),y(415:440),'-r','linewidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
grid on

% calculates the correlation coefficient b/w x and y in different time 
r = corrcoef(x,y)
r1 = corrcoef(x(377:388),y(377:388))
r2 = corrcoef(x(903:920),y(903:920))
r3 = corrcoef(x(415:440),y(415:440))
fileID = fopen('x.txt','w');
fprintf(fileID,'%f\n ',x);
fclose(fileID);
fileID = fopen('y.txt','w');
fprintf(fileID,'%f\n ',y);
fclose(fileID);

%A=load('x.txt');
% T=1:1:24001; %Time steps
% figure(2);
% plot(n,x,'-b','linewidth',2)
% hold on;
% plot(n,y,'-r','linewidth',2)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% grid on
% hold on;
% plot(T,A,'-black', 'markersize',2)% Plotting the time series for all realization