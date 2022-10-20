%% Q8 Part a
clearvars
clc

n_Vec = [200 2000];
K_Vec = [15 150];

for i=1:length(n_Vec)
    n = n_Vec(i);
    K = K_Vec(i);
    
    f = (1:n)';
    intcon = 1;

    A = [];
    b = [];

    Aeq = ones(n,1)';
    beq = K;

    lb = zeros(n,1);
    ub = ones(n,1)*2; 
    
    tic
    x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);
    toc
end
%% Q8 Part b
clearvars
clc

n_Vec = [200 2000];
K_Vec = [15 150];

for i=1:length(n_Vec)
    n = n_Vec(i);
    K = K_Vec(i);
    
    f = (1:n)';

    A = [];
    b = [];

    Aeq = ones(n,1)';
    beq = K;

    lb = zeros(n,1);
    ub = ones(n,1)*2; 
    
    tic
    x = linprog(f,A,b,Aeq,beq,lb,ub);
    toc
end