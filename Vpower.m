%% !Note: This is the earliest version of V-Power, just for reference.
%% Please DON'T run this script, run `Vpower2.m` instead.
% Virius Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power)
% author: Kun Yang

clc
clear
% [M]=xlsread('Mutation','MMFF');
% [P]=xlsread('Mutation','PPFF');
% [M] = readtable('MMFF.txt', 'Delimiter', 'tab');
% [P] = readtable('PPFF.txt', 'Delimiter', 'tab');
% M = M(:, 1:end);
% P = P(:, 1:end);
[M] = dlmread('MMFF_100.tsv', '\t', 1, 1);
[P] = dlmread('PPFF.txt', '\t', 1, 1);
M = M(:, 1:100);

SM = size(M);
SP = size(P);
X_all = zeros(SM(2), SP(2));
X0 = ones(SM(2), 1)./SM(2);
Aeq = ones(1, SM(2));
beq = 1;

LB = zeros(SM(2), 1);
UB = ones(SM(2), 1);
opts = optimoptions('fmincon','MaxFunctionEvaluations',300000,'MaxIterations',50000);
tic
for i = 1:SP(2)
    PP = P(:, i);
    [X, fval] = fmincon(@Preva,X0,[],[],Aeq,beq,LB,UB,[],opts,M,PP);
    X_all(:, i) = X;
end
toc

function f=Preva(X, MM, PP)
f=sum(sum(abs(MM * X - PP)));
end


