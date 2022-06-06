function [Population,FrontNo,fDS] = EnvironmentalSelection(Population,N,nb,mod)
% The environmental selection of DNEA

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Nojima, Y., Masuyama, N. and Shang, K., 2018, 
% September. A double-niched evolutionary algorithm and its behavior on 
% polygon-based problems. In International Conference on Parallel Problem 
% Solving from Nature (pp. 262-273). Springer, Cham.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    [Choose,fDS] = LastSelection(Population(Last).objs,Population(Last).decs,N-sum(Next),nb,mod);

    % Population for next generation
    fDS   = [zeros(1,sum(Next)),fDS(Choose)];
    Next(Last(Choose)) = true;
    Population = Population(Next);
    FrontNo    = FrontNo(Next);

end

function [Choose,fDS] = LastSelection(PopObj,PopDec,K,nb,mod)
% Select part of the solutions in the last front

    N = size(PopObj,1);    
    nb = min(nb,max(N-1,0));
    
    % Automatically calculate niche radius
    d_obj = pdist2(PopObj,PopObj,'euclidean');
    d_dec = pdist2(PopDec,PopDec,'euclidean');    
    temp = sort(d_obj);
    sigma_obj = sum(sum(temp(1:nb+1,:)))./(nb.*N);
    temp = sort(d_dec);
    sigma_dec = sum(sum(temp(1:nb+1,:)))./(nb.*N);
    
    if sigma_obj == 0
        sigma_obj = inf;
    end
    if sigma_dec == 0
        sigma_dec = inf;
    end
        
    neighbor_obj = d_obj < sigma_obj;
    neighbor_dec = d_dec < sigma_dec;
    
    Sh_obj = 1-d_obj./sigma_obj-eye(N);
    Sh_dec = 1-d_dec./sigma_dec-eye(N);
    
    % Calculate the double-sharing function of each solution
    switch mod
        case 0
            fDS = sum(Sh_obj.*neighbor_obj) + sum(Sh_dec.*neighbor_dec);
        case 1
            fDS = sum(Sh_obj.*neighbor_obj);
        case 2
            fDS = sum(Sh_dec.*neighbor_dec);
    end
            
    % Remove the solution with largest double-sharing function value
    Choose = true(1,N);
    while sum(Choose) > K
        [~,Del] = max(fDS); 
        Choose(Del) = false;
        neighbor_obj(Del,:) = false;
        neighbor_obj(:,Del) = false;
        neighbor_dec(Del,:) = false;
        neighbor_dec(:,Del) = false;
        switch mod
            case 0
                fDS = sum(Sh_obj.*neighbor_obj) + sum(Sh_dec.*neighbor_dec);
            case 1
                fDS = sum(Sh_obj.*neighbor_obj);
            case 2
                fDS = sum(Sh_dec.*neighbor_dec);
        end
        fDS(~Choose)=-inf;
    end
       
end