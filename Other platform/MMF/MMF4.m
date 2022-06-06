function varargout = MMF4(Operation,Global,input)
% <problem> <MMF>
% Multi-modal Multi-objective test Function
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright 2017-2018 Yiping Liu
% This is the code of MMF used in "Yiping Liu, Gary G. Yen, 
% and Dunwei Gong, A Multi-Modal Multi-Objective Evolutionary Algorithm 
% Using Two-Archive and Recombination Strategies, IEEE Transactions on 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% MMF is proposed in " Caitong Yue, Boyang Qu, and Jing Liang, 
% A Multi-objective Particle Swarm Optimizer Using Ring Topology for 
% Solving Multimodal Multi-objective Problems, IEEE Transactions on 
% Evolutionary Computation, 2017, Early Access".
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 2;
            Global.lower    = [-1 0];
            Global.upper    = [1 2];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {Global};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            PopDec(:,2) = PopDec(:,2)-(PopDec(:,2)>1);
            PopObj(:,1) = abs(PopDec(:,1));
            PopObj(:,2) = 1 - PopDec(:,1).^2 + 2.*(PopDec(:,2)-sin(pi.*PopObj(:,1))).^2;        
           
            PopCon = [];
           
            varargout = {PopObj};
        case 'PF'
            f =load('MMF4_Reference_PSPF_data');
            varargout = {f.PF};            
    end
end