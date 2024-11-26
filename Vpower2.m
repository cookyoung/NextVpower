% Virus Phylogenetic Resolver of Wastewater-based Epidemiology, version 2.0 (V_POWER2)
% author: Kun Yang


clc
clear

tic                                           % Record the time consumption of data reading
% [A,LA]=xlsread('MatrixA2024','A');                % Read the original "barcode" matrix A from the Excel file MatrixA
% [num1,SiteA]=xlsread('MatrixA2024','Sites');      % Read the mutation sites vector
% [P,SiteW]=xlsread('PFF2024','P');                % Read the mutation probability vectors P
% [num2,Samples]=xlsread('PFF2024','Samples');     % Read all sample names
A_table=dlmread('usher_barcodes.csv', ',', 0, 0);
A=A_table(2:end,2:end);
LA=A_table(2:end,1);
SiteA=A_table(1,2:end)';
P_table=dlmread('PP_raw_example.tsv', '\t', 0, 0);
P=P_table(2:end,2:end);
SiteW=P_table(2:end,1)';
Samples=P_table(1,2:end);
toc

% Filter mutation sites which were detected in wastewater but not recorded in matrix A
[isAW,posAW]=ismember(SiteW,SiteA);
K0=find(posAW==0);
% Record possible novel "key" mutation sites detected in sewage
NovelSites=SiteW(K0);  % Unrecorded mutation sites in matrix A
PNovel=P(K0,:);        % Probability of unrecorded mutation sites
% Pick up novel "key" mutation sites detected in most sewage samples (at least "2.0" eqivalent samples) but
% not yet recorded in matrix A
% "2.0" is an adjustable parameters indicating the mutation site detected
% in at list "2.0" eqivalent sewage samples
NovelK=find(sum(PNovel,2)>=2.0); 
NovelSitesK=NovelSites(NovelK);           % NovelSitesK record those novel "key" mutations are not yet recorded in the barcode matrix A
PNovelK=PNovel(NovelK,:);                 % PNovelK record the probability of novel "key" mutations detected in sewage samples.
%----------------------------------------------------------------
SiteW(K0)=[];               
posAW(K0)=[];
P(K0,:)=[];


% Filter some “old” lineages with fewer than 20 mutation sites 
% "20" is an adjustable parameters   
SUMAR=sum(A,2);
M=find(SUMAR<=20);
A(M,:)=[];
LA(M)=[];


% Filter lineages with mutation sites undetected in sewage samples
SUMAC=sum(A);
% Retain "key" mutation sites present in more than 200 lineages,
% "200" is an adjustable parameter
MuKey=find(SUMAC>=200);
PosAW=[posAW;MuKey'];
PosAW=sort(PosAW);
PosAW=unique(PosAW);
SSA=(1:size(SiteA));
SSA(:,PosAW)=[];
tic
for i=1:size(SSA,2)
    k=find(A(:,SSA(i))==1);
    A(k,:)=[];
    LA(k)=[];
end
toc

% Filter mutation sites undetected in sewage sample
ssA=(1:size(SiteA));
ssA(:,posAW)=[];
A(:,ssA)=[];
SiteA(ssA)=[];


% Filter mutation sites not recorded in "Barcode" matrix A 
% Repeated process,can be omitted,just for verification
[isWA,posWA]=ismember(SiteA,SiteW);
SSW=(1:size(SiteW));
SSW(posWA)=[];
P(SSW,:)=[];
SiteW(SSW)=[];


% Deconvolution
SA = size(A);
SP = size(P);
X_all = zeros(SA(1), SP(2));
X0 = ones(SA(1), 1)./SA(1);
Aeq = ones(1, SA(1));
beq = 1;
LB = zeros(SA(1), 1);
UB = ones(SA(1), 1);

tic
for i = 1:SP(2)
    PP = P(:, i);
    [X, fval] = fmincon(@Preva,X0,[],[],Aeq,beq,LB,UB,[],[],A',PP);
    X_all(:, i) = X;
end
toc

AA=bar(X_all',0.4,'stacked');

function f=Preva(X, MM, PP)
f=sum(sum(abs(MM * X - PP)));
end









    








