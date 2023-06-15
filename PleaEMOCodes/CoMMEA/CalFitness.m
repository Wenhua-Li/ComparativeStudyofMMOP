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

if strcmp(Operation,'LocalC') % Dividing all solution into niches
    mod = 1;
    dist = pdist2(PopDec,PopDec);
    if mod==0
        R = mean(mean(dist))/2;
        C=[];
        for i=1:N
            C(i).p = i; C(i).decs = Population(i).decs';
        end
        while 1
            decs = [C.decs]';
            dist = pdist2(decs,decs);
            dist(logical(eye(length(C))))=inf;
            [a,b] = find(dist==min(min(dist)));
            p1=a(1);p2=a(2);
            if dist(p1,p2)>R
                break;
            end
            C(p1).p=[C(p1).p C(p2).p];
            C(p1).decs = mean(Population(C(p1).p).decs)';
            C(p2)=[];
        end
       
        K=length(C);
        Niching=false(N);
        for i=1:K
            cluster=C(i).p;
            for j=1:length(cluster)
                Niching(cluster,j) = true;
                Niching(j,cluster) = true;
            end
        end
        fprintf('n:%.2f N:%d AvgN: %.1f \n',R,N,mean(sum(Niching)))
    elseif mod == 1
        %% Method 1: distance based
        if size(PopDec,2)<=8
            R = mean(mean(dist))/4;
        else
            R = mean(mean(dist))/2;
        end
%         R = .4*prod(max(PopDec)-min(PopDec)).^(1./size(PopDec,2));
%         R = 0.5*prod(max(PopDec)-min(PopDec)).^(1./size(PopDec,2));
        Niching = dist<R;
%         fprintf('R:%.2f N:%d AvgN: %.1f \n',R,N,mean(sum(Niching)))
    end
    Dominate = Dominate & Niching;
end

%% Calculate S(i)
S = sum(Dominate,2);

%% Calculate R(i)
R = zeros(1,N);
for i = 1 : N
    R(i) = sum(S(Dominate(:,i)));
end

%% Calculate D(i)
% Distance = pdist2(PopDec,PopDec);
% Distance(logical(eye(length(Distance)))) = inf;
% Distance = sort(Distance,2);
% D = 1./(Distance(:,floor(sqrt(N)))+2);

%% Calculate the fitnesses
% Fitness = R + D';
Fitness = R;
end