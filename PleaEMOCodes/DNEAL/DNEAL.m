function DNEAL(Global)
% <algorithm> <A>
% Searching for Local Pareto Optimal Solutions: A Case Study on
% Polygon-Based Problems
% mod --- 0 --- 0 both, 1 obj, 2 dec
% nb --- 3 --- Neighbor size
% K --- 4 --- the number of non-dominated fronts to store
% Nns --- 50 --- the number of nearest solutions to be considered

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Nojima, Y., Masuyama, N. and Han, Y., 2019, June.
% Searching for local pareto optimal solutions: A case study on
% polygon-based problems. In 2019 IEEE Congress on Evolutionary Computation
% (CEC) (pp. 896-903). IEEE.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

%% Parameter setting
[mod,nb,K,Nns] =  Global.ParameterSet(0,3,3,50);

%% Generate random population
Population = Global.Initialization();
[~,FrontNo,fDS] = MultiFront(Population,Global.N,nb,mod,K,Nns);

%% Optimization
while Global.NotTermination(Population)
    numPop = length(Population);
    FrontNo(FrontNo<=K)=1;
    MatingPool = TournamentSelection(2,numPop,FrontNo,fDS);
    Offspring  = Global.Variation(Population(MatingPool(1:Global.N)));
    [Population,FrontNo,fDS] = MultiFront([Population,Offspring],Global.N,nb,mod,K,Nns);
end
end
