function Fitness = CalFitness(PopObj,PopDec,state)
% Calculate the fitness of each solution

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

N = size(PopObj,1);

%% Detect the dominance relation between each two solutions
Dominate = false(N);
for i = 1 : N-1
    for j = i+1 : N
        k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
        if k == 1
            Dominate(i,j) = true;
        elseif k == -1
            Dominate(j,i) = true;
        end
    end
end

%     [dist,index] = sort(pdist2(PopDec,PopDec),2);
dist=pdist2(PopDec,PopDec);
V=.5*prod(max(PopDec)-min(PopDec)).^(1./size(PopDec,2));
S = zeros(N,1);
for i=1:N
    index = (dist(i,:)<V);
    if sum(index)<0
        S(i) = sum(Dominate(i,:));
    else
        S(i) = sum(Dominate(i,index));
    end
end

%     n = round(sqrt(N/2)/2)
%     Dominates = Dominate()
%% Calculate S(i)
%     nDominate = Dominate(index);nDominate=nDominate(:,1:n);
%     S = sum(nDominate,2);

%% Calculate R(i)

nn = min(N/40,state*N/20);

R = zeros(1,N);
for i = 1 : N
    index = (dist(i,:)<V);
    if sum(index)<nn
        R(i) = sum(S(Dominate(:,i)));
    else
        R(i) = sum(S(Dominate(index,i)));
    end
end
%% Calculate D(i)
Distance = pdist2(PopDec,PopDec);
Distance(logical(eye(length(Distance)))) = inf;
Distance = sort(Distance,2);
D = 1./(Distance(:,floor(sqrt(N)))+2);

%% Calculate the fitnesses
Fitness = R + D';
end