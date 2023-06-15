function varargout = PolygonProblem(Operation,Global,input)
% <problem> <PolygonProblem>
% The multi-modal multi-objective polygon distance minimization problem
%
% lower --- -50 --- Lower bound of decision variables
% upper ---  50 --- Upper bound of decision variables
% row   ---   2   --- Number of polygons in a row
% col   ---  2   --- Number of polygons in a column
% distance ---  5 --- Distance between polygons center
[lower, upper, obj.row, obj.col, obj.distance] = Global.ParameterSet(-50, 50, 2, 2, 5);
Polygons = CreatePolygons(obj.row, obj.col, obj.distance, Global.M);
switch Operation
    case 'init'
        Global.M        = 6;
        Global.D        = 4;
        Global.lower    = lower+zeros(1,Global.D);
        Global.upper    = upper+ones(1,Global.D);
        Global.operator = @EAreal;
        PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
        %rand(input,Global.D);
        varargout = {PopDec};
    case 'value'
        N = obj.row * obj.col;
        PopDec = input;
        % Min distance to point i in each polygon
        PopObj = Inf;
        for i=1:N
            % Transform 2d to D-dimensional plane using two basis (0, 1, 0, 1, ...) and (1, 0, 1, 0, ...).
            vertices = repmat(Polygons(:, :, i), 1, Global.D / 2);
            
            val = pdist2(PopDec, vertices);
            PopObj = min(PopObj, val);
        end
        
        PopCon = [];
        
        varargout = {input,PopObj,PopCon};
    case 'PF'
        [PS, PF] = PSnPF(Polygons, Global.D, Global.N);
%         f = UniformPoint(input,Global.M);
        varargout = {PS, PF};
end
end