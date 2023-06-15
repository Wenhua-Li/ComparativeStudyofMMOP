function [Population,Fitness] = EnvironmentalSelection1(Population,N,Operation,EvoState)
% The environmental selection of CoMMEA

%% SPEA2 selection
Fitness = CalFitness(Population,Operation);
Next = Fitness < 1;
if sum(Next) < N
    [~,Rank] = sortrows([Fitness' -Crowding(Population.decs)]);
    Population = Population(Rank(1:N));
else
    Population = Population(Next);
    while length(Population)>N
        dist = sort(pdist2(Population.decs,Population.decs));
        CrowdDis = sum(dist(1:3,:));
        [~,index] = min(CrowdDis);
        Population(index) = [];
    end
end
if EvoState<0.5
    Fitness=CalFitness(Population,Operation)-Crowding(Population.decs)';
else
    Fitness=-Crowding(Population.decs)';
end

% %% IBEA selection: Better performance on many-objective problems
% kappa = 0.05;
% [Fitness,I,C] = CalFitness3(Population,kappa);
% Next = 1 : length(Population);
% while length(Next) > N
%     [~,x]   = min(Fitness(Next));
%     Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/kappa);
%     Next(x) = [];
% end
% Population = Population(Next);
% Fitness=Kdis(Population,N/2);
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