function varargout = MMF12(Operation,Global,input)
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
            Global.lower    = [0 0];
            Global.upper    = [1 1];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {Global};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            PopObj(:,1) = PopDec(:,1);
            q=4;
            Alfa=2;
            num_of_peak=2;
            g=2-((sin(num_of_peak*pi.*PopDec(:,2))).^6).*(exp(-2*log10(2).*((PopDec(:,2)-0.1)/0.8).^2));
            h=1-(PopObj(:,1)./g).^Alfa-(PopObj(:,1)./g).*sin(2*pi*q*PopObj(:,1));
            PopObj(:,2) = g.*h;
            
            PopCon = [];
    
            varargout = {PopObj};
        case 'PF'
            f =load('MMF12_Reference_PSPF_data');
            varargout = {f.PF};        
    end
end