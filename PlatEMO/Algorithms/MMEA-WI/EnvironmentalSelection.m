function [Population,Fitness] = EnvironmentalSelection(Population,N,kappa,state)
% The environmental selection of IBEA

Next = 1 : length(Population);
[Fitness,I,C] = CalFitness(Population,kappa);

dist = pdist2(Population.decs,Population.decs);

var=Population.decs;

sigma=prod(max(var)-min(var))^(1/size(var,2))/size(var,1);
sigma=sigma*(exp(-state));
% weight = 1./(1+exp(dist./sigma));

weight=1/(sigma*(2*pi)^0.5);
weight=weight.*exp(-dist.^2/(2*sigma^2));

newfit=zeros(size(Fitness));
for i=1:size(var,1)
    tmp=weight(i,:).*(Fitness);
    newfit(i)=sum(tmp);
end

weight=sum(weight);
newfit=newfit.*weight;

% fd=zeros(size(var,1));
% for i=1:size(var,1)
%     for j=1:size(var,1)
%         fd(i,j)=dist(i,j)/(1+0.5*(newfit(i)+newfit(j)));
%     end
% end
% Fitness=1./(1+sum(fd));

Fitness = newfit;

while length(Next) > N
    [~,x]   = min(Fitness(Next));
    Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/kappa);
    Next(x) = [];
end
Population = Population(Next);
% Fitness=Fitness(Next);

Fitness=Kdis(Population,N/2);

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
end