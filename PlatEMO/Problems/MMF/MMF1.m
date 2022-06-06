function varargout = MMF1(Operation,Global,input)
% <problem> <MMF>
% Multi-modal Multi-objective test Function
% operator --- EAreal

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 2;
            Global.lower    = [1 -1];
            Global.upper    = [3 1];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            PopObj(:,1) = abs(PopDec(:,1)-2);
            PopObj(:,2) = 1-sqrt(PopObj(:,1))+2*(PopDec(:,2)-sin(6*pi*PopObj(:,1)+pi)).^2;
            
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f =load('MMF1_Reference_PSPF_data');
            varargout = {f.PF};        
    end
end