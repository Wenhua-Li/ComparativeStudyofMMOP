function [Population,Fitness] = EnvironmentalSelection2(Population,N,Operation,eps)
% The environmental selection of CoMMEA

%% Epsilon-dominance-based method
[FrontNo,~] = NDSort(Population.objs,inf);
Next = find(FrontNo == 1);
Remain = find(FrontNo > 1);
isEpsSol = false(1,numel(Remain));

for j=1:numel(Remain)
    TmpObj = [Population(Remain(j)).objs; (1+eps).*Population(Next).objs];
    [fno,~] = NDSort(TmpObj,inf);
    if fno(1)==1
        isEpsSol(j) = true;
    end
end

if length(Next) + sum(isEpsSol) < N
%     fprintf('Way1 PF:%3d  EpsSol:%3d \n',length(Next), sum(isEpsSol))
    Cdis = -Crowding(Population.decs);
    [~,idx] = sortrows([FrontNo' Cdis],'ascend');
    Population = Population(idx(1:N));
else
%     fprintf('Way2 PF:%3d  EpsSol:%3d \n',length(Next), sum(isEpsSol))
    Population = [Population(Next) Population(Remain(isEpsSol))];
end

%% Calculate local convergence quality
Fitness = CalFitness(Population,Operation);

%% Crowding distance based second selection
Next = Fitness < 1;
K = 3;
% fprintf('LocalPF:%3d \n',sum(Next))
if sum(Next) < N
    dist = sort(pdist2(Population.decs,Population.decs));
    CrowdDis = -sum(dist(1:K,:));
    [~,idx] = sortrows([Fitness' CrowdDis']);
    Population = Population(idx(1:N));
else
    Population = Population(Next);
    %% way1
    while length(Population)>N
        dist = sort(pdist2(Population.decs,Population.decs));
        CrowdDis = sum(dist(1:K,:));
        [~,index] = min(CrowdDis);
        Population(index) = [];
    end
    
%     %% way2
%     while length(Population)>N
%         dobj = pdist2(Population.objs,Population.objs); dobj = dobj./max(max(dobj));
%         ddec = pdist2(Population.decs,Population.decs); ddec = ddec./max(max(ddec));
%         dist = sort(dobj)/2+sort(ddec);
%         CrowdDis = sum(dist(1:K,:));
%         [~,index] = min(CrowdDis);
%         Population(index) = [];
%     end
end

%% Calculate fitness for parents selection
% Fitness=-Crowding(Population.decs)-Crowding(Population.objs);
% Fitness=-Crowding(Population.decs);
% dist = sort(pdist2(Population.decs,Population.decs));
% CrowdDis = sum(dist(1:K,:));
% Fitness = -CrowdDis;
Fitness = Kdis(Population,N/2);
end

function fDN = Kdis(Pop,K)
PopObj=Pop.objs;
PopDec=Pop.decs;
Np = size(PopObj,1);
d_dec = pdist2(PopDec,PopDec,'euclidean');
d_dec(logical(eye(Np))) = inf;
sdd = sort(d_dec);
dn_dec = sum(sdd(1:K,:));
avg_dn_dec = mean(dn_dec);
if avg_dn_dec == 0
    avg_dn_dec = inf;
end
fDN = 1./(1+dn_dec./avg_dn_dec);
end