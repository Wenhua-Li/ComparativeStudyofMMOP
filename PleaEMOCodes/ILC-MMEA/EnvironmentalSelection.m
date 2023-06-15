function [Population,Fitness] = EnvironmentalSelection(Population,N)
% The environmental selection of ILC_MMEA

%% Calculate the fitness of each solution
Fitness = CalFitness(Population.objs,Population.decs);

%% Environmental selection
Next = Fitness < 1;
if sum(Next) < N
    [~,Rank] = sort(Fitness);
    Next(Rank(1:N)) = true;
elseif sum(Next) > N
    Del  = Truncation(Population(Next).decs,Population(Next).objs,sum(Next)-N);
    Temp = find(Next);
    Next(Temp(Del)) = false;
end
% Population for next generation
Population = Population(Next);
Fitness    = Fitness(Next);
end

function Del = Truncation(PopDec,PopObj,K)
% Select part of the solutions by truncation
%% Truncation
Del = false(1,size(PopObj,1));
while sum(Del) < K
    Remain   = find(~Del);
    Distance = Crowding(PopDec(Remain,:));%+Crowding(PopObj(Remain,:));
    [~,Rank] = min(Distance);
    Del(Remain(Rank(1))) = true;
end
end

function [CrowdDis]=Crowding(Pop)
%%Harmonic average distance of each solution in the decision space
% return: the crowding distance of each individual
[N, ~]=size(Pop);
K=N-1;
Z = min(Pop,[],1);
Zmax = max(Pop,[],1);
pop=(Pop-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
distance=pdist2(pop,pop);
[value,~]=sort(distance,2);
CrowdDis=K./sum(1./value(:,2:N),2);
end