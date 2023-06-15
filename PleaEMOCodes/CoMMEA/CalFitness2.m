function Fitness = CalFitness(Population,Operation)
% Calculate the fitness of each solution

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

PopObj = Population.objs;
PopDec = Population.decs;
N = size(PopObj,1);

dist = pdist2(PopDec,PopDec);
R = mean(mean(dist))/4;
Niching = dist<R;
rank = zeros(1,N);
for i=1:N
    p = find(Niching(i,:)==1);
    [FrontNo,~] = NDSort(Population(p).objs,length(p));
    rank(p) = rank(p)+FrontNo./length(p);    
end

Fitness = rank;
end