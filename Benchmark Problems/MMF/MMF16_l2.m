function varargout = MMF16_l2(Operation,Global,input)
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
            Global.M        = 3;
            Global.D        = 3;
            Global.lower    = zeros(1,Global.D);
            Global.upper    = ones(1,Global.D);
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            M=Global.M;
            num_of_peak=2;
            x=PopDec;
            
            
                num_of_g_peak=1;%   number of global PSs
    num_of_l_peak=2;%   number of local PSs


% if size(x,1)~=1
%     x=reshape(x,1,length(x));
% end
%  
for i=1:size(x,1)
   
    if x(i,end)>=0&&x(i,end)<0.5
        g(i,:)=2-(sin(2*num_of_g_peak*pi.*x(i,end))).^2;%
    elseif x(i,end)>=0.5&&x(i,end)<=1
        g(i,:)=2-exp(-2*log10(2).*((x(i,end)-0.1)/0.8).^2).*(sin(2*num_of_l_peak*pi.*x(i,end))).^2;
    end
end
 y = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,M-1:-1:1)*pi/2)];
        
%             
%             g=2-exp(-2*log10(2).*((x(:,end)-0.1)/0.8).^2).*(sin(num_of_peak*pi.*x(:,end))).^2;
%             y = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(x(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(x(:,M-1:-1:1)*pi/2)];
            PopObj = y;
            
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f =load('MMF15_Reference_PSPF_data');
            varargout = {f.PF};        
    end
end