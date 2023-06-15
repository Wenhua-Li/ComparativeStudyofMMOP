function [ps,pf]=SSMOPSO(func_name,VRmin,VRmax,n_obj,Particle_Number,Max_Gen,var)


% Dimension: n_var --- dimensions of decision space
%            n_obj --- dimensions of objective space
%% Input:
%                      Dimension                    Description
%      func_name       1 x length function name     the name of test function     
%      VRmin           1 x n_var                    low bound of decision variable
%      VRmax           1 x n_var                    up bound of decision variable
%      n_obj           1 x 1                        dimensions of objective space
%      Particle_Number 1 x 1                        population size
%      Max_Gen         1 x 1                        maximum  generations

%% Output:
%                     Description
%      ps             Pareto set
%      pf             Pareto front


%% Initialize parameters
    n_var=size(VRmin,2);              
    Max_FES=Max_Gen*Particle_Number;   
    n_PBA=40;     
    Rs=0.1;
    cc=[2.05 2.05];                    
    iwt=0.7298;
%% Initialize particles' positions and velocities
    mv=0.5*(VRmax-VRmin);
    VRmin=repmat(VRmin,Particle_Number,1);
    VRmax=repmat(VRmax,Particle_Number,1);
    Vmin=repmat(-mv,Particle_Number,1);
    Vmax=-Vmin;
    pos=VRmin+(VRmax-VRmin).*rand(Particle_Number,n_var); 
    vel=Vmin+2.*Vmax.*rand(Particle_Number,n_var);        
%% Evaluate the population
    fitness=zeros(Particle_Number,n_obj);
    for i=1:Particle_Number;
        fitness(i,:)=feval(func_name,'value',var,pos(i,:));
    end
    fitcount=Particle_Number;            
    particle=[pos,fitness];              
    row_of_cell=ones(1,Particle_Number); 
    col_of_cell=size(particle,2);       
    PBA=mat2cell(particle,row_of_cell,col_of_cell);
    p=[];
for i=1:Max_Gen
    fitpop=non_domination_scd_sort(particle(:,1:n_var+n_obj), n_obj, n_var);
    spop=[];
    for a=1:Particle_Number
        for b=1:Particle_Number
            if particle(a,1:n_var)== fitpop(b,1:n_var);
                p(a)=fitpop(b,n_var+n_obj+1);
                fitpop(b,1:n_var)=1000.*ones(1,n_var);%[1000 1000];  
                break;  
            end
        end
    end
      [sortp,sortindex]=sort(p,'ascend');
     popsort=particle(sortindex,1:n_var);
     popvel=vel(sortindex,1:n_var);
     psort=p(sortindex);
     u=1;
   while u<=600&size(popsort,1)~=0 
    dist=zeros(size(popsort,1),1);
    dist(1:size(popsort,1),:)=sum((ones(size(popsort,1),1)*popsort(1,:)-popsort).^2,2)<Rs.^2;
    spop(u).pop=popsort(dist==1,:);
    spop(u).vel=popvel(dist==1,:);
   spop(u).species=popsort(1,:);
    popsort=popsort(dist==0,:);
     popvel= popvel(dist==0,:);
    u=u+1;
   end
    for t=1:size(spop,2)      
         for g=1:size(spop(t).pop,1)
             for h=1:Particle_Number
                 if spop(t).pop(g,:)==particle(h,1:n_var)
                     v=h;
                     particle(h,1:n_var)=1000.*ones(1,n_var);%[1000 1000];
                     break;
                 end
             end
             PBA_v=PBA{v,1};      
            pbest=PBA_v(1,:);    
                 spop(t).vel(g,1:n_var)=iwt.*spop(t).vel(g,1:n_var)+cc(1).*rand(1,n_var).*(pbest(1,1:n_var)-spop(t).pop(g,1:n_var))+cc(2).*rand(1,n_var).*( spop(t).species(1,1:n_var)-spop(t).pop(g,1:n_var));
                spop(t).vel(g,1:n_var)=(spop(t).vel(g,1:n_var)>mv).*mv+(spop(t).vel(g,1:n_var)<=mv).*spop(t).vel(g,1:n_var); 
           spop(t).vel(g,1:n_var)=(spop(t).vel(g,1:n_var)<(-mv)).*(-mv)+(spop(t).vel(g,1:n_var)>=(-mv)).*spop(t).vel(g,1:n_var);
              spop(t).pop(g,1:n_var)=spop(t).vel(g,1:n_var)+spop(t).pop(g,1:n_var);
                  
         spop(t).pop(g,1:n_var)=((  spop(t).pop(g,1:n_var)>=VRmin(1,:))&(  spop(t).pop(g,1:n_var)<=VRmax(1,:))).*  spop(t).pop(g,1:n_var)...
               +(  spop(t).pop(g,1:n_var)<VRmin(1,:)).*(VRmin(1,:)+0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,n_var))+(  spop(t).pop(g,1:n_var)>VRmax(1,:)).*(VRmax(1,:)-0.25.*(VRmax(1,:)-VRmin(1,:)).*rand(1,n_var));
                 
              particle(v,1:n_var)=spop(t).pop(g,:);
              vel(v,1:n_var)=spop(t).vel(g,1:n_var);
              
            fitness(v,:)=feval(func_name,'value',var,spop(t).pop(g,:));
            fitcount=fitcount+1;
            particle(v,1:n_var+n_obj)=[spop(t).pop(g,:),fitness(v,:)];
            PBA_v=[PBA_v(:,1:n_var+n_obj);particle(v,:)]; 
             PBA_v = non_domination_scd_sort(PBA_v(:,1:n_var+n_obj), n_obj, n_var);
             if size(PBA_v,1)>n_PBA
                 PBA{v,1}=PBA_v(1:n_PBA,:);
             else
                 PBA{v,1}=PBA_v;
             end
              if fitcount>Max_FES
                break;
              end
  
         end
          if fitcount>Max_FES
                break;
           end
    end
      
   
   if fitcount>Max_FES
        break;
    end
   
end


%% Output ps and pf
    tempEXA=cell2mat(PBA); 
    tempEXA=non_domination_scd_sort(tempEXA(:,1:n_var+n_obj), n_obj, n_var);
     if size(tempEXA,1)>Particle_Number
         EXA=tempEXA(1:Particle_Number,:);
     else
        EXA=tempEXA;
     end
   tempindex=find(EXA(:,n_var+n_obj+1)==1);
   ps=EXA(tempindex,1:n_var);
   pf=EXA(tempindex,n_var+1:n_var+n_obj);
end



function f = non_domination_scd_sort(x, n_obj, n_var)

%% Input£º
%                      Dimension                      Description
%      x               num_particle x n_var+n_obj     population to be sorted     
%      n_obj           1 x 1                          dimensions of objective space
%      n_var           1 x 1                          dimensions of decision space

%% Output:
%              Dimension                                  Description
%      f       N_particle x (n_var+n_obj+4)               Sorted population  
%    in f      the (n_var+n_obj+1)_th column stores the front number
%              the (n_var+n_obj+2)_th column stores the special crowding distance   
%              the (n_var+n_obj+3)_th column stores the crowding distance in decision space
%              the (n_var+n_obj+4)_th column stores the crowding distance in objective space


    [N_particle, ~] = size(x);% Obtain the number of particles

% Initialize the front number to 1.
    front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
    F(front).f = [];
    individual = [];

%% Non-Dominated sort. 

    for i = 1 : N_particle
        % Number of individuals that dominate this individual
        individual(i).n = 0; 
        % Individuals which this individual dominate
        individual(i).p = [];
        for j = 1 : N_particle
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k = 1 : n_obj
                if (x(i,n_var + k) < x(j,n_var + k))
                    dom_less = dom_less + 1;
                elseif (x(i,n_var + k) == x(j,n_var + k))  
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if dom_less == 0 && dom_equal ~= n_obj
                individual(i).n = individual(i).n + 1;
            elseif dom_more == 0 && dom_equal ~= n_obj
                individual(i).p = [individual(i).p j];
            end
        end   
        if individual(i).n == 0
            x(i,n_obj + n_var + 1) = 1;
            F(front).f = [F(front).f i];
        end
    end
% Find the subsequent fronts
    while ~isempty(F(front).f)
       Q = [];
       for i = 1 : length(F(front).f)
           if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ...
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        x(individual(F(front).f(i)).p(j),n_obj + n_var + 1) = ...
                            front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
               end
           end
       end
       front =  front + 1;
       F(front).f = Q;
    end
% Sort the population according to the front number
    [~,index_of_fronts] = sort(x(:,n_obj + n_var + 1));
    for i = 1 : length(index_of_fronts)
        sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
    end
    current_index = 0;

%% SCD. Special Crowding Distance

    for front = 1 : (length(F) - 1)
  
        crowd_dist_obj = 0;
        y = [];
        previous_index = current_index + 1;
        for i = 1 : length(F(front).f)
            y(i,:) = sorted_based_on_front(current_index + i,:);
        end
        current_index = current_index + i;
   % Sort each individual based on the objective
        sorted_based_on_objective = [];
        for i = 1 : n_obj+n_var
            [sorted_based_on_objective, index_of_objectives] = ...
                sort(y(:,i));
            sorted_based_on_objective = [];
            for j = 1 : length(index_of_objectives)
                sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            end
            f_max = ...
                sorted_based_on_objective(length(index_of_objectives), i);
            f_min = sorted_based_on_objective(1,  i);

            if length(index_of_objectives)==1
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;  
            elseif i>n_var
              
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;
                y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)=0;
            else
              
                 y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)...
                    = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
                sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min);
                 y(index_of_objectives(1),n_obj + n_var + 1 + i)=2*(sorted_based_on_objective(2, i)-...
                sorted_based_on_objective(1, i))/(f_max - f_min);
            end
             for j = 2 : length(index_of_objectives) - 1
                next_obj  = sorted_based_on_objective(j + 1, i);
                previous_obj  = sorted_based_on_objective(j - 1,i);
                if (f_max - f_min == 0)
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = 1;
                else
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = ...
                         (next_obj - previous_obj)/(f_max - f_min);
                end
             end
        end
    %% Calculate distance in decision space
        crowd_dist_var = [];
        crowd_dist_var(:,1) = zeros(length(F(front).f),1);
        for i = 1 : n_var
            crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,n_obj + n_var + 1 + i);
        end
        crowd_dist_var=crowd_dist_var./n_var;
        avg_crowd_dist_var=mean(crowd_dist_var);
    %% Calculate distance in objective space
        crowd_dist_obj = [];
        crowd_dist_obj(:,1) = zeros(length(F(front).f),1);
        for i = 1 : n_obj
            crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,n_obj + n_var + 1+n_var + i);
        end
        crowd_dist_obj=crowd_dist_obj./n_obj;
        avg_crowd_dist_obj=mean(crowd_dist_obj);
    %% Calculate special crowding distance
        special_crowd_dist=zeros(length(F(front).f),1);
        for i = 1 : length(F(front).f)
            if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
                special_crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i)); 
            else
                special_crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i)); 
            end
        end
        y(:,n_obj + n_var + 2) = special_crowd_dist;
        y(:,n_obj+n_var+3)=crowd_dist_var;
        y(:,n_obj+n_var+4)=crowd_dist_obj;
        [~,index_sorted_based_crowddist]=sort(special_crowd_dist,'descend');
        y=y(index_sorted_based_crowddist,:);
        y = y(:,1 : n_obj + n_var+4 );
        z(previous_index:current_index,:) = y;
    end
    
f = z();
end


