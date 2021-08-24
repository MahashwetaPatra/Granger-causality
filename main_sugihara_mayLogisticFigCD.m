% NOTE: This code performs Sugihara-May 1990 algorithm, using an ensemble
%       mean method for prediction model instead of the weighted sum as in 
%       the original method.
%
% HIST: - Jan 11, 2021: Created by Weiran Cai from CK's algorithm
%       - Feb 02, 2021: Fixed a bug in index searching by CK
%       - Feb 11, 2021: added more instructions by CK 
%       - Jun 29, 2021: added weighed sum method by Patra
%       - Jul 01, 2021: fixed bugs in sum/weight calculation by CK 
%       - Aug 10,2021: shows that cross mapping of Y using M_X fails; however, the cross map of X succeeds by Patra
%==========================================================================
%
% Generate logistic output data arrays now
%
clc;close all;clear all;
x(1) =0.4;
y(1) =0.2;
n=1:1:24001;
for j = 1:24000
      x(j+1) = x(j)*(3.8-3.8*x(j)-0.0*y(j));
      y(j+1) = y(j)*(3.5-3.5*y(j)-0.1*x(j));
end
z=x;
x=y;
y=z;
%x=20*sin(2*pi/200*(1:24001));
%x=load('x.txt');
%y=load('y.txt'); 
%
% reading an input data or geneating it from a given function
% %
l_array=[];
rho_array=[];
l_train=2000;
%for l_train = 50:50:4000 
    %
    % define a training set and all required parameters
    % 
    A_train = x(1:l_train);           % training data set
    %l_train = length(A_train);       % length of training set
    m = 3;                            % embedding dim
    tau = 1;                         % delay time, previously it was tau=10 
    num_nb = 4;                       % number of neighbours
    X_em = cstr_em(A_train, m, tau);  % construct the phase space
    Tp = 1;                           % the maximum forecast lead time
    l_test = 1500;                    % length of test set, previoisly it was 2000
    A_fwd = zeros(l_test,Tp);         % array of prediction for test set
    %
    % generating a prediction for all test data in the test set A\A_train.
    %
    t_pr0 = 100;
    for t_pr = (t_pr0+1):(t_pr0+l_test)
        %
        % searching for the indices of all neighbors in the training set for 
        % each predictee
        %
        x_pr = zeros(1,3);
        x_pr(1,:) = x((t_pr-(m-1)*tau):tau:t_pr);
        [ID_ss] = find_ss(x_pr, X_em, m, l_train, num_nb,tau);
        %
        % produce prediction for the entire range of lead time from 1-Tp
        %
        x1=[x(ID_ss(1)), x(ID_ss(1)-tau), x(ID_ss(1)-2*tau)];
        x2=[x(ID_ss(2)), x(ID_ss(2)-tau), x(ID_ss(2)-2*tau)];
        x3=[x(ID_ss(3)), x(ID_ss(3)-tau), x(ID_ss(3)-2*tau)];
        x4=[x(ID_ss(4)), x(ID_ss(4)-tau), x(ID_ss(4)-2*tau)];
        u1=exp(-norm(x_pr(1,:)-x1)/norm(x_pr(1,:)-x1));
        u2=exp(-norm(x_pr(1,:)-x2)/norm(x_pr(1,:)-x1));
        u3=exp(-norm(x_pr(1,:)-x3)/norm(x_pr(1,:)-x1));
        u4=exp(-norm(x_pr(1,:)-x4)/norm(x_pr(1,:)-x1));              
        sum=u1+u2+u3+u4;
        w(1)=u1/sum;
        w(2)=u2/sum;
        w(3)=u3/sum;
        w(4)=u4/sum;
        %w(1)+w(2)+w(3)+w(4);
        %A_fwd(t_pr-t_pr0,t_fwd) = mean(A(ID_ss+t_fwd));
        t_fwd=1;
        sum1=0.0;
        for i=1:4            
            %sum1=sum1+w(i)*x(ID_ss(i)+t_fwd);
            sum1=sum1+w(i)*y(ID_ss(i)+t_fwd);
            %sum1=sum1+0.25*A(ID_ss(i)+t_fwd);
        end
        A_fwd(t_pr-t_pr0,t_fwd) = sum1;        
    end
    %
    % compute correlation between forecast and actual value for a single lead 
    % times t_fwd = 1
    %
    %rho_corr = zeros(Tp,1);
    t_fwd=1;    
    R = (corrcoef(A_fwd(:,t_fwd),y((1:l_test)+t_pr0+t_fwd)));
    %R = (corrcoef(A_fwd(:,t_fwd),x((1:l_test)+t_pr0+t_fwd)));
    rho_corr = R(1,2)
    l_array=[l_array;l_train];
    rho_array=[rho_array;rho_corr];
%end
figure(1);
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontSize',14,'LineWidth',2.75)
plot(l_array,rho_array,'.-r', 'markersize',20)% Plotting the time series for all realization

figure(2);
set(gca, 'GridLineStyle', ':') %dotted grid lines
set(gca,'FontSize',14,'LineWidth',2.75)
plot(A_fwd(:,t_fwd),y((1:l_test)+t_pr0+t_fwd),'+r', 'markersize',10)% Plotting the time series for all realization
xlim([0 1])
ylim([0 1])
