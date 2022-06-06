function MMEAWI(Global)
% <algorithm> <A>
% Indicator-Based Selection in Multiobjective Search
% kappa --- 0.05 --- Fitness scaling factor

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Parameter setting
kappa = Global.ParameterSet(0.05);
t_gen=ceil(Global.maxgen*0.4);

%% Generate random population
Population = Global.Initialization();
[Population,pfit] = EnvironmentalSelection(Population,Global.N,kappa,1);
[Arcd,afit] = UpdateArc(Population,Population,Global.N);

%% Optimization
while Global.NotTermination(Arcd)
    if Global.gen>=t_gen && rand>0.5 % stage 2
        [~,p1]=min(afit);
        joint=[Arcd Population];
        dist=pdist2(Arcd(p1).decs,joint.decs);
        [~,so]=sort(dist,'ascend');
        MatingPool= so(randperm(round(Global.N/5),round(Global.N/10))+1);
        parents=[Arcd(p1) joint(MatingPool)];
        Offspring  = Global.Variation(parents);
        [Population,pfit] = EnvironmentalSelection([Population,Offspring],Global.N,kappa,Global.gen./Global.maxgen);
        [Arcd,afit] = UpdateArc(Arcd,Offspring,Global.N);
    else %stage 1
        MatingPool = TournamentSelection(2,round(Global.N/10),pfit);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,pfit] = EnvironmentalSelection([Population,Offspring],Global.N,kappa,Global.gen./Global.maxgen);
        [Arcd,afit] = UpdateArc(Arcd,Offspring,Global.N);
    end
%     if Global.gen==Global.maxgen
%        save('r.mat','Population','Arcd');
%     end
end
end