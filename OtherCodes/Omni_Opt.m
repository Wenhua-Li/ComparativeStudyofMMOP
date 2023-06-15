function [ps,pf]=Omni_Opt(func_name,xl,xu,n_obj,pop,Maxgen,var)
n=size(xl,2);               %Obtain the dimensions of decision space
feval_max=pop*Maxgen;
feval_count=0;

for i=1:pop
    particle(i,1:n) = xl+(xu-xl).*rand(1,n);
end

for i=1:pop
%     particle(i,n+1:n+n_obj)=feval(func_name,particle(i,1:n));
    particle(i,n+1:n+n_obj)=feval(func_name,'value',var,particle(i,1:n));
    feval_count=feval_count+1;%%每个个体评价一次，feval_count就加1
end

particle(:,1:n+n_obj+2) = non_domination_sort_crowd_dist(particle, n_obj, n);

for i=1:Maxgen
    
%    disp(['Omni_opt  gen=  ' num2str(i) ]);
    
    pool = round(pop/2);
    tour = 2;
    
    parent_chromosome = tournament_selection(particle, pool, tour);
    
    mu = 20;
    mum = 20;
    offspring_chromosome = genetic_operator(parent_chromosome, n_obj, n, mu, mum, xl, xu,func_name,feval_count,var);
    
    [main_pop,temp] = size(particle);
    [offspring_pop,temp] = size(offspring_chromosome);
    clear temp
    
    intermediate_chromosome(1:main_pop,:) = particle;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : n+n_obj) = offspring_chromosome;
    
    intermediate_chromosome = non_domination_sort_crowd_dist(intermediate_chromosome, n_obj, n);
    
    particle = replace_chromosome(intermediate_chromosome, n_obj, n, pop);
    
    if feval_count>feval_max,break, end

end
ps=particle(:,1:n);
pf=particle(:,n+1:n+n_obj);

end

function f = non_domination_sort_crowd_dist(x, M, V)

%% function f = non_domination_sort_crowd_dist(x, M, V)
% x     粒子
% M     目标空间维度
% N     决策空间维度

% f     按照非支配排序从前到后，拥挤距离从大到小排序后的种群
%%
[N, m] = size(x);
clear m

% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = [];

%% Non-Dominated sort. 

for i = 1 : N
    % Number of individuals that dominate this individual
    individual(i).n = 0; 
    % Individuals which this individual dominate
    individual(i).p = [];
    for j = 1 : N
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : M
            if (x(i,V + k) < x(j,V + k))
                dom_less = dom_less + 1;
            elseif (x(i,V + k) == x(j,V + k))  
                dom_equal = dom_equal + 1;
            else
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= M
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= M
            individual(i).p = [individual(i).p j];
        end
    end   
    if individual(i).n == 0
        x(i,M + V + 1) = 1;
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
               		x(individual(F(front).f(i)).p(j),M + V + 1) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end

[temp,index_of_fronts] = sort(x(:,M + V + 1));
for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
end
current_index = 0;

%% Crowding distance

for front = 1 : (length(F) - 1)
%    objective = [];
    crowd_dist_obj = 0;
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);%取出第front前沿个体放入y中
    end
    current_index = current_index + i;
    % Sort each individual based on the objective
    sorted_based_on_objective = [];
    for i = 1 : M+V
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
            y(index_of_objectives(1),M + V + 1 + i) = 1;%如果该前沿只有一个点
          
        elseif i>V
              %目标空间内每一维的边界点拥挤距离赋值为最大值或最小值（最小化问题低边界点赋最大值高边界点赋最小值，最大化问题颠倒）
            y(index_of_objectives(1),M + V + 1 + i) = 1;
            y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)=0;
        else
            %决策空间内每一维的边界点拥挤距离赋值为：边界点与邻居点距离的二倍
             y(index_of_objectives(length(index_of_objectives)),M + V + 1 + i)...
                = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
            sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min);
             y(index_of_objectives(1),M + V + 1 + i)=2*(sorted_based_on_objective(2, i)-...
            sorted_based_on_objective(1, i))/(f_max - f_min);
        end
         for j = 2 : length(index_of_objectives) - 1
            next_obj  = sorted_based_on_objective(j + 1, i);
            previous_obj  = sorted_based_on_objective(j - 1,i);
            if (f_max - f_min == 0)
                y(index_of_objectives(j),M + V + 1 + i) = 1;
            else
                y(index_of_objectives(j),M + V + 1 + i) = ...
                     (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    %Calculate distance in x space
    crowd_dist_var = [];
    crowd_dist_var(:,1) = zeros(length(F(front).f),1);
    for i = 1 : V
        crowd_dist_var(:,1) = crowd_dist_var(:,1) + y(:,M + V + 1 + i);
    end
    crowd_dist_var=crowd_dist_var./V;
    avg_crowd_dist_var=mean(crowd_dist_var);
    %Calculate distance in f space
    crowd_dist_obj = [];
    crowd_dist_obj(:,1) = zeros(length(F(front).f),1);
    for i = 1 : M
        crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,M + V + 1+V + i);
    end
    crowd_dist_obj=crowd_dist_obj./M;
    avg_crowd_dist_obj=mean(crowd_dist_obj);
    crowd_dist=zeros(length(F(front).f),1);
    for i = 1 : length(F(front).f)
        if crowd_dist_obj(i)>avg_crowd_dist_obj||crowd_dist_var(i)>avg_crowd_dist_var
            crowd_dist(i)=max(crowd_dist_obj(i),crowd_dist_var(i));
        else
            crowd_dist(i)=min(crowd_dist_obj(i),crowd_dist_var(i));
        end
    end
    y(:,M + V + 2) = crowd_dist;
    y(:,M+V+3)=crowd_dist_var;
    y(:,M+V+4)=crowd_dist_obj;
    [~,index_sorted_based_crowddist]=sort(crowd_dist,'descend');%%同一非支配排序的再按照拥挤距离排序
    y=y(index_sorted_based_crowddist,:);
    y = y(:,1 : M + V+2 );
    z(previous_index:current_index,:) = y;
end
f = z();
end
function f = tournament_selection(chromosome, pool_size, tour_size)


[pop, variables] = size(chromosome);
% The peunltimate element contains the information about rank.
rank = variables - 1;
% The last element contains information about crowding distance.
distance = variables;

% Until the mating pool is filled, perform tournament selection
for i = 1 : pool_size
    % Select n individuals at random, where n = tour_size
    for j = 1 : tour_size
        % Select an individual at random
        candidate(j) = round(pop*rand(1));
        % Make sure that the array starts from one. 
        if candidate(j) == 0
            candidate(j) = 1;
        end
        if j > 1
            % Make sure that same candidate is not choosen.
            while ~isempty(find(candidate(1 : j - 1) == candidate(j)))
                candidate(j) = round(pop*rand(1));
                if candidate(j) == 0
                    candidate(j) = 1;
                end
            end
        end
    end
    % Collect information about the selected candidates.
    for j = 1 : tour_size
        c_obj_rank(j) = chromosome(candidate(j),rank);
        c_obj_distance(j) = chromosome(candidate(j),distance);
    end
    % Find the candidate with the least rank
    min_candidate = ...
        find(c_obj_rank == min(c_obj_rank));
    % If more than one candiate have the least rank then find the candidate
    % within that group having the maximum crowding distance.
    if length(min_candidate) ~= 1
        max_candidate = ...
        find(c_obj_distance(min_candidate) == max(c_obj_distance(min_candidate)));
        % If a few individuals have the least rank and have maximum crowding
        % distance, select only one individual (not at random). 
        if length(max_candidate) ~= 1
            max_candidate = max_candidate(1);
        end
        % Add the selected individual to the mating pool
        f(i,:) = chromosome(candidate(min_candidate(max_candidate)),:);
    else
        % Add the selected individual to the mating pool
        f(i,:) = chromosome(candidate(min_candidate(1)),:);
    end
end
end
function f  = replace_chromosome(intermediate_chromosome, M, V,pop)

[N, m] = size(intermediate_chromosome);

% Get the index for the population sort based on the rank
[temp,index] = sort(intermediate_chromosome(:,M + V + 1));

clear temp m

% Now sort the individuals based on the index
for i = 1 : N
    sorted_chromosome(i,:) = intermediate_chromosome(index(i),:);
end

% Find the maximum rank in the current population
max_rank = max(intermediate_chromosome(:,M + V + 1));

% Start adding each front based on rank and crowing distance until the
% whole population is filled.
previous_index = 0;
for i = 1 : max_rank
    % Get the index for current rank i.e the last the last element in the
    % sorted_chromosome with rank i. 
    current_index = max(find(sorted_chromosome(:,M + V + 1) == i));
    % Check to see if the population is filled if all the individuals with
    % rank i is added to the population. 
    if current_index > pop
        % If so then find the number of individuals with in with current
        % rank i.
        remaining = pop - previous_index;
        % Get information about the individuals in the current rank i.
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        % Sort the individuals with rank i in the descending order based on
        % the crowding distance.
        [temp_sort,temp_sort_index] = ...
            sort(temp_pop(:, M + V + 2),'descend');
        % Start filling individuals into the population in descending order
        % until the population is filled.
        for j = 1 : remaining
            f(previous_index + j,:) = temp_pop(temp_sort_index(j),:);
        end
        return;
    elseif current_index < pop
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
    else
        % Add all the individuals with rank i into the population.
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end
    % Get the index for the last added individual.
    previous_index = current_index;
end
end
function f  = genetic_operator(parent_chromosome, M, V, mu, mum, l_limit, u_limit,func_name,feval_count,var)


[N,m] = size(parent_chromosome);

clear m
p = 1;
% Flags used to set if crossover and mutation were actually performed. 
was_crossover = 0;
was_mutation = 0;


for i = 1 : N
    % With 90 % probability perform crossover
    if rand(1) < 0.9
        % Initialize the children to be null vector.
        child_1 = [];
        child_2 = [];
        % Select the first parent
        parent_1 = round(N*rand(1));
        if parent_1 < 1
            parent_1 = 1;
        end
        % Select the second parent
        parent_2 = round(N*rand(1));
        if parent_2 < 1
            parent_2 = 1;
        end
        % Make sure both the parents are not the same. 
        while isequal(parent_chromosome(parent_1,:),parent_chromosome(parent_2,:))
            parent_2 = round(N*rand(1));
            if parent_2 < 1
                parent_2 = 1;
            end
        end
        % Get the chromosome information for each randomnly selected
        % parents
        parent_1 = parent_chromosome(parent_1,:);
        parent_2 = parent_chromosome(parent_2,:);
        % Perform corssover for each decision variable in the chromosome.
        for j = 1 : V
            % SBX (Simulated Binary Crossover).
            % For more information about SBX refer the enclosed pdf file.
            % Generate a random number
            u(j) = rand(1);
            if u(j) <= 0.5
                bq(j) = (2*u(j))^(1/(mu+1));
            else
                bq(j) = (1/(2*(1 - u(j))))^(1/(mu+1));
            end
            % Generate the jth element of first child
            child_1(j) = ...
                0.5*(((1 + bq(j))*parent_1(j)) + (1 - bq(j))*parent_2(j));
            % Generate the jth element of second child
            child_2(j) = ...
                0.5*(((1 - bq(j))*parent_1(j)) + (1 + bq(j))*parent_2(j));
            % Make sure that the generated element is within the specified
            % decision space else set it to the appropriate extrema.
            if child_1(j) > u_limit(j)
                child_1(j) = u_limit(j);
            elseif child_1(j) < l_limit(j)
                child_1(j) = l_limit(j);
            end
            if child_2(j) > u_limit(j)
                child_2(j) = u_limit(j);
            elseif child_2(j) < l_limit(j)
                child_2(j) = l_limit(j);
            end
        end
        % Evaluate the objective function for the offsprings and as before
        % concatenate the offspring chromosome with objective value.
%         child_1(:,V + 1: M + V) = evaluate_objective(child_1, M, V);
%         child_2(:,V + 1: M + V) = evaluate_objective(child_2, M, V);
          child_1(:,V + 1: M + V) = feval(func_name,'value',var,child_1(1:V));
          child_2(:,V + 1: M + V) = feval(func_name,'value',var,child_2(1:V));%feval(func_name,child_2(1:V));
          
          feval_count=feval_count+2;


        % Set the crossover flag. When crossover is performed two children
        % are generate, while when mutation is performed only only child is
        % generated.
        was_crossover = 1;
        was_mutation = 0;
    % With 10 % probability perform mutation. Mutation is based on
    % polynomial mutation. 
    else
        % Select at random the parent.
        parent_3 = round(N*rand(1));
        if parent_3 < 1
            parent_3 = 1;
        end
        % Get the chromosome information for the randomnly selected parent.
        child_3 = parent_chromosome(parent_3,:);
        % Perform mutation on eact element of the selected parent.
        for j = 1 : V
           r(j) = rand(1);
           if r(j) < 0.5
               delta(j) = (2*r(j))^(1/(mum+1)) - 1;
           else
               delta(j) = 1 - (2*(1 - r(j)))^(1/(mum+1));
           end
           % Generate the corresponding child element.
           child_3(j) = child_3(j) + delta(j);
           % Make sure that the generated element is within the decision
           % space.
           if child_3(j) > u_limit(j)
               child_3(j) = u_limit(j);
           elseif child_3(j) < l_limit(j)
               child_3(j) = l_limit(j);
           end
        end
        % Evaluate the objective function for the offspring and as before
        % concatenate the offspring chromosome with objective value.    
%         child_3(:,V + 1: M + V) = evaluate_objective(child_3, M, V);
        child_3(:,V + 1: M + V) = feval(func_name,'value',var,child_3(1:V));%feval(func_name,child_3(1:V));
        feval_count=feval_count+1;
        % Set the mutation flag
        was_mutation = 1;
        was_crossover = 0;
    end
    % Keep proper count and appropriately fill the child variable with all
    % the generated children for the particular generation.
    if was_crossover
        child(p,:) = child_1;
        child(p+1,:) = child_2;
        was_cossover = 0;
        p = p + 2;
    elseif was_mutation
        child(p,:) = child_3(1,1 : M + V);
        was_mutation = 0;
        p = p + 1;
    end
end
f = child;
end




