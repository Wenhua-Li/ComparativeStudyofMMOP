function HREA(Global)
% <algorithm> <A>
% Hierarchy Ranking Based Multimodal Multi-objective Evolutionary Algorithm
% eps --- 0.3 --- parameter for quality of the local PF
% p --- 0.5 --- parameter for probability

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
[eps,p] = Global.ParameterSet(0.3,0.5);

%% Generate random population
Population = Global.Initialization();
[~,CrowdDis1] = MPFEnvironmentalSelection(Population,Global.N);
[Archive,CrowdDis2] = ArchiveUpdate(Population,Global.N,eps,0);

%% Optimization
while Global.NotTermination(Archive)
    if Global.evaluated >= Global.evaluation * 0.5 && rand < p
        MatingPool2 = TournamentSelection(2,round(Global.N),-CrowdDis2);
        Offspring  = Global.Variation([Archive(MatingPool2)]);
    else
        MatingPool1 = TournamentSelection(2,round(Global.N),-CrowdDis1);
        Offspring  = Global.Variation([Population(MatingPool1)]);
    end
    [Population,CrowdDis1] = MPFEnvironmentalSelection([Population,Offspring],Global.N);
    [Archive,CrowdDis2] = ArchiveUpdate([Archive,Offspring],Global.N,eps,Global.evaluated / Global.evaluation);

    save('HREA_mmf15a.mat','Archive')
end

end