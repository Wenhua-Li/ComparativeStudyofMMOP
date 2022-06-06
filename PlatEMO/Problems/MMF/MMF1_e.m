function varargout = MMF1_e(Operation,Global,input)
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
            Global.lower    = [1 -20];
            Global.upper    = [3 20];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            PopObj(:,1) = abs(PopDec(:,1)-2);
            index1=PopDec(:,1)<2;
            PopObj(index1,2) = 1-sqrt(PopObj(index1,1))+2*(PopDec(index1,2)-sin(6*pi*PopObj(index1,1)+pi)).^2;
            index2=PopDec(:,1)>=2;
            PopObj(index2,2) = 1-sqrt(PopObj(index2,1))+2*(PopDec(index2,2)-exp(PopDec(index2,1)).*sin(6*pi*PopObj(index2,1)+pi)).^2;
            
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f =load('MMF1_e_Reference_PSPF_data');
            varargout = {f.PF};        
    end
end