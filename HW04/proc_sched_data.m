% Data for minimum energy processor speed scheduling.
clearvars
clc

n = 12;  % number of jobs.
T = 16;  % number of time periods.

Smin = 1;  % min processor speed.
Smax = 4;  % max processor speed.
R = 1;  % max slew rate.

% Parameters in power/speed model.
alpha = 1;
beta = 1;
gamma = 1;

% Job arrival times and deadlines.
A = [1; 3;  4; 6; 9;  9; 11; 12; 13; 13; 14; 15];
D = [3; 6; 11; 7; 9; 12; 13; 15; 15; 16; 14; 16];
% Total work for each job.
W = [2; 4; 10; 2; 3; 2; 3; 2; 3; 4; 1; 4];

% Plot showing job availability times & deadlines.
figure;
hold on;
scatter(A,[1:n],'k*');
scatter(D+1,[1:n],'ko');
for i=1:n
    plot([A(i),D(i)+1],[i,i],'k-');
end
hold off;
xlabel('time t');
ylabel('job i');

%%
% th n T
% S T
cvx_begin
    variables X(n,T) S(T)
    minimize( phi(S, alpha, beta, gamma) )
    subject to
        for i=1:T
            S(i) <= Smax
            S(i) >= Smin
        end
        for i=1:T-1
            S(i) - S(i+1) <= R
            S(i+1) - S(i) <= R
        end
        for i=1:n
            sum( X(i,A(i):D(i)) ) >= W(i)
        end
        for i=1:T
            ones(1,n) * X(:,i) == S(i)
        end
        X > 0
cvx_end

for i=1:T
    X(:,i) = X(:,i)/S(i);
end
%% Plot showing job availability times & deadlines & u efforts
addpath(genpath('matlabGiftiCifti'))
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(3,1,1)
hold on;
scatter(A,[1:n],'k*');
scatter(D+1,[1:n],'ko');
for i=1:n
    plot([A(i),D(i)+1],[i,i],'black');
end
for i=1:n
    [x,y] = find(X(i,:) > 0.005);
    plot([y(1),y(end)+1],[i,i],'r');
end
hold off;
xlabel('Time');
ylabel('Job');
L = {'Job Available' 'Working on the Job'};
xlim([1, 16])
legend(L,'Location','bestoutside')
subplot(3,1,2)
plot(X')
box off
xlabel('Time')
ylabel('\theta')
L = {};
for i=1:n
    L{i} = "Job " + num2str(i);
end
legend(L,'Location','bestoutside')
xlim([1, 16])
subplot(3,1,3)
bar((S*ones(1,n)).*X',1,'stacked')
L = {};
for i=1:n
    L{i} = "Job " + num2str(i);
end
legend(L,'Location','bestoutside')
export_fig('Q1.png','-r600');

phi(S, alpha, beta, gamma)