function f = non_domination_scd_kmeans_sort(x, n_obj, n_var,KK)
% non_domination_scd_sort:  sort the population according to non-dominated relationship and special crowding distance

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
%% 
% perfrom k-means in each front, then calculate the decision space distance
finalpop = [];
avg_crowd_dist_var = [];
for front = 1 : (length(F)-1)
    front_index = find(sorted_based_on_front(:,n_obj + n_var + 1)==front);
    zhongjianpop = sorted_based_on_front(front_index,:);
    K = ceil(length(front_index)/KK);
    if size(zhongjianpop,1) < K
        keyboard
    end
    [ind,~] = kmeans(zhongjianpop(:,1:n_var),K);
    pop=cell(max(ind),1);
    for i1=1:length(ind)
        pop{ind(i1)} = [pop{ind(i1)};zhongjianpop(i1,:)];
    end
    
    for i2= 1 : size(pop,1)
        crowd_dist_obj = 0;
        y = pop{i2};
        sorted_based_on_objective = [];
        for i = 1 : n_var
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
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;  %If there is only one point in current front
            elseif length(index_of_objectives)==2
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;
                y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)=1;
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
        crowd_dist_var(:,1) = zeros(size(pop{i2},1),1);
        for ii = 1 : n_var
            crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,n_obj + n_var + 1 + ii);
        end
        crowd_dist_var=crowd_dist_var./n_var;
        %avg_crowd_dist_var_k=mean(crowd_dist_var);
        y(:,n_obj+n_var+2)=crowd_dist_var;
        y = y(:,1 : n_obj + n_var+2 );
        finalpop=[finalpop;y];
       % avg_crowd_dist_var(front) = [avg_crowd_dist_var(front);avg_crowd_dist_var_k];
    end
end

sorted_based_on_front = [];
%% calculate the objective space distance
[~,index_of_fronts] = sort(finalpop(:,n_obj + n_var + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = finalpop(index_of_fronts(i),:);
   % avg_crowd_dist_var(i,:) = avg_crowd_dist_var(index_of_fronts(i),:);
end

  current_index = 0;

%% CSCD

    for front = 1 : (length(F) - 1)
  
        crowd_dist_obj = 0;
        y = [];
        previous_index = current_index + 1;
        for i = 1 : length(F(front).f)
            y(i,:) = sorted_based_on_front(current_index + i,:);%put the front_th rank into y
        end
        current_index = current_index + i;
   % Sort each individual based on the objective
        sorted_based_on_objective = [];
        for i = n_var+1 : n_obj+n_var
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
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;  %If there is only one point in current front
            elseif i>n_var
                % deal with boundary points in objective space
                % In minimization problem, set the largest distance to the low boundary points and the smallest distance to the up boundary points
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;
                y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)=0;
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
        avg_crowd_dist_var = mean(y(:,n_var+n_obj+2));
        for i = 1 : length(F(front).f)
            crowd_dist_var(i) = y(i,n_var+n_obj+2);
            if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
                special_crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i)); % Eq. (6) in the paper
            else
                special_crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i)); % Eq. (7) in the paper
            end
        end
        y(:,n_obj + n_var + 3) = special_crowd_dist;
        y(:,n_obj+n_var+4)=crowd_dist_obj;
        [~,index_sorted_based_crowddist]=sort(special_crowd_dist,'descend');%sort the particles in the same front according to CSCD
        y=y(index_sorted_based_crowddist,:);
        y = y(:,1 : n_obj + n_var+4 );
        z(previous_index:current_index,:) = y;
    end
f= z();
end
