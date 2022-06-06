function MMOEAC(Global)
% <algorithm> <A>
% date: 2019.3.22(final version)
% author: Wu Lin
% Please contact {linwu2017@email.szu.edu.cn} if you have any problem.

%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    delta=5;
    
    %% Generate random population
    Population = Global.Initialization();
    gen=1;
    %% Optimization
    while Global.NotTermination(Population)
        % MatingPool:tournament selection based on decision spatial density
        CrowdDis=Crowding(Population.decs);
        MatingPool=TournamentSelection(2,Global.N,-CrowdDis);
        % Recombination
        Offspring  = Global.Variation(Population(MatingPool));
        Union=[Population Offspring];
        % Environmental Selection
        Population=Environmental_Selection(Union,Global.N,delta);
        gen=gen+1;
        
%         if mod(Global.gen,5)==0
%             save(['Data\mmoeac_' num2str(Global.gen) '.mat'],'Population');
%         end
        
    end
end

