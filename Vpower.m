% Virius Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power)
% (Pipeline Version)
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
M = M(:, 1:89);
%MFF=M(1:25,:);
SM = size(M);
SP = size(P);
X_all = zeros(SM(2), SP(2));
X0 = ones(SM(2), 1)./SM(2);
Aeq = ones(1, SM(2));
beq = 1;
% A=Aeq;
% b=beq;
LB = zeros(SM(2), 1);
UB = ones(SM(2), 1);
% opts = optimoptions('fmincon','MaxFunctionEvaluations',300000,'MaxIterations',50000);
tic
for i = 1:SP(2)
    PP = P(:, i);
%     [X, fval] = fmincon(@Preva,X0,[],[],Aeq,beq,LB,UB,[],[],M,PP,opts);
%     [X, fval] = fmincon(@Preva,X0,M,PP,Aeq,beq,LB,UB,[],opts);
    [X, fval] = fmincon(@Preva,X0,[],[],Aeq,beq,LB,UB,[],[],M,PP);
    X_all(:, i) = X;
end
%XX=MFF'\P;
% X0=ones(SM(1),SP(2))./SM(1);
%X01=ones(1,SP(2));
%X02=zeros(SM(1)-1,SP(2));
%X0=[X01;X02];
%[X,fval]=fminsearch(@Preva,X0,[],MFF,P);
toc

function f=Preva(X, MM, PP)
f=sum(sum(abs(MM * X - PP)));
end


