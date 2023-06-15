function CoMMEA2(Global)
% <algorithm> <A>
% CoMMEA: Coevolutionary Framework for Generalized Multimodal Multi-objective Optimization
% eps --- 0.3 --- parameter for controlling the quality of the local PF

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

eps = Global.ParameterSet(0.3);
%% Generate random population
Population1 = Global.Initialization();
Population2 = Global.Initialization();
Fitness1    = CalFitness(Population1,'Normal');
Fitness2    = CalFitness(Population2,'LocalC');

%% Optimization
while Global.NotTermination(Population2)
    MatingPool1 = TournamentSelection(2,round(Global.N/2),Fitness1);
    MatingPool2 = TournamentSelection(2,Global.N,Fitness2);
    Offspring1  = Global.Variation(Population1(MatingPool1));
    Offspring2  = Global.Variation(Population2(MatingPool2));
    
    %% Enviornmental selection
    EvoState = Global.evaluated/Global.evaluation;
    CurEps = max(-2*log2(2*EvoState),eps); % Compute the current epsilon value
    [Population1,Fitness1] = EnvironmentalSelection1([Population1,Offspring1,Offspring2],Global.N,'Normal',EvoState);
    [Population2,Fitness2] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],Global.N,'LocalC',CurEps);
end

end

%         hold on
%         PopObj = Population1.objs;
%         PopDec = Population1.decs;
%         if size(PopObj,2)==3
%             plot3(PopObj(:,1),PopObj(:,2),PopObj(:,3),'r.')
%         elseif size(PopObj,2)==2
%             plot(PopObj(:,1),PopObj(:,2),'r.')
%         else
%             plot(PopObj,'r')
%         end