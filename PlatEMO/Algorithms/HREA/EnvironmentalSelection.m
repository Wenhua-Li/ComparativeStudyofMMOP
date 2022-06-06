function [Population,CrowdDis] = EnvironmentalSelection(Population,N)

n = length(Population);

%% Cal the local convergence
dist = pdist2(Population.decs,Population.decs);
V=0.2*prod(max(Population.decs)-min(Population.decs)).^(1./size(Population.decs,2));

DominationX = zeros(n); % Pareto domination relationship between pair of solutions
for i=1:n
    for j=i+1:n
        if dist(i,j) > V
            continue
        end
        L1 = Population(i).objs < Population(j).objs;
        L2 = Population(i).objs > Population(j).objs;
        if all(L1|(~L2))
            DominationX(i,j) = 0;
            DominationX(j,i) = 1;
        elseif all(L2|(~L1))
            DominationX(i,j) = 1;
            DominationX(j,i) = 0;
        end
    end
end

LocalC = zeros(1,n);
for i=1:n
    tmp = dist(i,:);
    index = tmp<V;
    LocalC(i) = (sum(DominationX(i,index)))./sum(index);
end
% CrowdDis = Crowding(Population.decs)'; % 用下面的拥挤距离效果更好一点

dist = sort(pdist2(Population.decs,Population.decs));
CrowdDis = sum(dist(1:3,:));

[~,index] = sortrows([LocalC' -CrowdDis']);
Population = Population(index);

if length(Population)>N
    Population = Population(1:N);
end

CrowdDis = Crowding(Population.decs);

end