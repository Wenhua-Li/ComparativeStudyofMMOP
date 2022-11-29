function f  = replace_decision_chromosome_kmeans(intermediate_chromosome, M, V,pop)

%第二判定标准选x空间距离，按x空间拥挤距离降序排列，选前几个
[ C,IA,IC ] = unique(intermediate_chromosome(:,1:V),'rows');
intermediate_chromosome = intermediate_chromosome(IA,:);
%           IA为矩阵C中的元素在矩阵A中的位置，
%
%           IC为矩阵A中的元素在矩阵C中的位置。

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
            sort(temp_pop(:, M + V + 3),'descend');
        %%
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


