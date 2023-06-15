function [Population,dk] = UpdateArc(Population,offspring,N)
% Update the Archive
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang,
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary
% Multi-Objective Optimization [Educational Forum], IEEE Computational
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

joint=[Population offspring];
% [Fitness,I,C] = CalFitness(joint,0.05);
% Next = 1 : length(joint);
% while length(Next) > N+5
%     [~,x]   = min(Fitness(Next));
%     Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/0.05);
%     Next(x) = [];
% end
% Population = joint(Next);
% length(Population)

[FrontNo,MaxFNo] = NDSort(joint.objs,N);
next=FrontNo==1;
Population=joint(next);

[Choose,fDKN] = DoubleNearestSelection(Population,N,5);
Population = Population(Choose);
dk = fDKN(Choose);

    function [Choose,fDN] = DoubleNearestSelection(Pop,N,K)
        PopObj=Pop.objs;
        PopDec=Pop.decs;
        Np = size(PopObj,1);
        Choose = true(1,Np);
        if Np <= K
            fDN = zeros(1,Np);
            return;
        end
        d_obj = pdist2(PopObj,PopObj,'euclidean');
        d_dec = pdist2(PopDec,PopDec,'euclidean');
        d_obj(logical(eye(Np))) = inf;
        d_dec(logical(eye(Np))) = inf;
        
        sdo = sort(d_obj);
        sdd = sort(d_dec);
        dn_obj = sum(sdo(1:K,:));
        dn_dec = sum(sdd(1:K,:));
        avg_dn_obj = mean(dn_obj);
        avg_dn_dec = mean(dn_dec);
        if avg_dn_obj == 0
            avg_dn_obj = inf;
        end
        if avg_dn_dec == 0
            avg_dn_dec = inf;
        end
        fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
        
        while sum(Choose) > N
            [~,Del] = max(fDN);
            
            Choose(Del) = false;
            d_obj(Del,:) = inf;
            d_obj(:,Del) = inf;
            d_obj(Del,:) = inf;
            d_obj(:,Del) = inf;
            
            sdo = sort(d_obj);
            sdd = sort(d_dec);
            dn_obj = sum(sdo(1:K,:));
            dn_dec = sum(sdd(1:K,:));
            avg_dn_obj = mean(dn_obj);
            avg_dn_dec = mean(dn_dec);
            if avg_dn_obj == 0
                avg_dn_obj = inf;
            end
            if avg_dn_dec == 0
                avg_dn_dec = inf;
            end
            fDN = 1./(1+dn_obj./avg_dn_obj+dn_dec./avg_dn_dec);
            fDN(~Choose) = -inf;
        end
    end
end
