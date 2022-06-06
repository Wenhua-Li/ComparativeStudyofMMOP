function [Population,FrontNo,fDS] = MultiFront(Population,N,nb,mod,K,Nns)
% The environmental selection of DNEAL

%------------------------------- Reference --------------------------------
% Liu, Y., Ishibuchi, H., Nojima, Y., Masuyama, N. and Han, Y., 2019, June.
% Searching for local pareto optimal solutions: A case study on 
% polygon-based problems. In 2019 IEEE Congress on Evolutionary Computation
% (CEC) (pp. 896-903). IEEE.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
    
    %% find local non-dominated solutions
    NP = length(Population);
    NB = min(Nns,max(NP-1,0));
    d_dec = pdist2(Population.decs,Population.decs,'euclidean');    
    temp = sort(d_dec);
    th_dec = sum(sum(temp(1:NB+1,:)))./(NB.*NP);       
    local = true(1,NP);
    for i = 1:NP
        Distance = d_dec(i,:);
        Nbs = find(Distance < th_dec);
        I = Nbs == i;
        [FNo,~] = NDSort(Population(Nbs).objs,inf);
        if(FNo(I)>1)
          local(i)=false;
        end
    end
       
    %% preserve the first K non-dominated fronts
    numL = sum(local);
    LPop = Population(local);
    [FrontNo,MaxFN] = NDSort(LPop.objs,inf); 
    Choose = FrontNo <= K;
    NumMF = sum(Choose);
    
    %% remove solutions with bad double-sharing function values
    if NumMF<N
        NLPop = Population(~local);
        FrontNo2 = NDSort(NLPop.objs,inf);
        FrontNo = [FrontNo,FrontNo2+N];
        i = 0;
        temp = 0;
        IFN = unique(FrontNo);
        while temp < N
            i = i+1;
            temp = temp + sum(FrontNo==IFN(i));
        end  
        MaxFNo = IFN(i);        
        Population = [LPop,NLPop];
        Next = FrontNo < MaxFNo;
        Last   = find(FrontNo==MaxFNo);
        [Choose,fDS] = DoulbeSharingSelection(Population(Last).objs,Population(Last).decs,N-sum(Next),nb,mod);        
        % Population for next generation
        fDS   = [zeros(1,sum(Next)),fDS(Choose)];
        Next(Last(Choose)) = true;
        Population = Population(Next);
        FrontNo    = FrontNo(Next);
        
    else
        fDS = zeros(1,numL);
        Choose = false(1,numL);
        for i=1:min(K,MaxFN)
            C1 = find(FrontNo==i);           
            Ni = length(C1);
            Ni = min(Ni,N);
            [C2,ds] = DoulbeSharingSelection(LPop(C1).objs,LPop(C1).decs,Ni,nb,mod);
            fDS(C1) = ds;
            C3=C1(C2);
            Choose(C3) = true;
        end
         % Population for next generation
        Population = LPop(Choose);
        fDS = fDS(Choose);
        FrontNo = FrontNo(Choose);
    end
    
    
   

end

function [Choose,fDS] = DoulbeSharingSelection(PopObj,PopDec,Ni,nb,mod)
% Select part of the solutions using the doulbe sharing function

    N = size(PopObj,1);    
    nb = min(nb,max(N-1,0));
    
    d_obj = pdist2(PopObj,PopObj,'euclidean');
    d_dec = pdist2(PopDec,PopDec,'euclidean');
    
    temp = sort(d_obj);
    sigma_obj = sum(sum(temp(1:nb+1,:)))./(nb.*N);
    temp = sort(d_dec);
    sigma_dec = sum(sum(temp(1:nb+1,:)))./(nb.*N);
    
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
               
    Choose = true(1,N);
    while sum(Choose) > Ni
        [~,Del] = max(fDS+Choose-1); 
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
    end
    
    
end