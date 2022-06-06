function Population=Environmental_Selection(Union, N, delta)
% The environmental selection
    
    %% Neighborhood-based Clustering Method
    range=max(Union.decs,[],1)-min(Union.decs,[],1);
    r=(range)*0.1;
    
    C=NCM(Union.decs,r);
    K=length(C);
    local_A=[];
    % local_A: reserve the nondominated solutions in each local clusters
    for i=1:K
        cluster=C(i).p;
        if length(cluster)<=delta
           continue;
        end
        [FrontNo1,~]=NDSort(Union(cluster).objs,length(cluster));
        local_A=[local_A cluster(FrontNo1==1)];
    end
    
    %% P: the number of individuals is not less than N
    [FrontNo,~]=NDSort(Union.objs,length(Union));
    P = find(FrontNo==1);
    temp=setdiff(P,local_A);
    P=[local_A temp];
    i=1;
    while  length(P)<N
        i=i+1;
        temp1=find(FrontNo==i);
        temp2=setdiff(temp1, local_A);
        P=[P temp2]; 
    end

    if length(P)>N 
        newpop=Union(P);
        z=min(newpop.objs,[],1);
        Z=max(newpop.objs,[],1);
        % normalization in the decision space
        ZZZ=(newpop.objs-z)./repmat(Z-z,length(newpop),1);
        % hierarchical clustering method
        H=clusterdata(ZZZ,'maxclust',N,'distance','euclidean','linkage','ward');
        
        count=length(newpop)-N;
        %% delete redundant individual in turn
        for i=1:count
            CrowdDis=Crowding(newpop.decs);
            num=hist(H,1:N);
            % find the most crowding clusters in the objective space 
            % and then add them into R
            I=find(num==max(num));
            R=[];
            for j=1:length(I)
                R=[R;find(H==I(j))];
            end
            % delete the solution with minimal decision spatial crowding distance 
            T=find(min(CrowdDis(R))==CrowdDis(R));
            s=randperm(length(T));
            x=R(T(s(1)));
            newpop(x)=[];
            H(x)= [];
        end
        Population=newpop;
    else
        Population=Union(P);
    end
end