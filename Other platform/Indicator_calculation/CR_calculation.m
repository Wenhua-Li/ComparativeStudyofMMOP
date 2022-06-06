function CR=CR_calculation(obtained_ps,reference_ps)
% CR_calculation: Calculate the cover rate of the obtained Pareto set

% Dimension: n_var --- dimensions of decision space

%% Input:
%                      Dimension                                                   Description
%      obtained_ps     population_size x n_var                                     Obtained Pareto set     
%      reference_ps    num_of_solutions_in_reference_ps x n_var                    Reference Pareto set

%% Output:
%                     Description
%      CR             Cover rate of the obtained Pareto set

n_var=size(reference_ps,2);
%find the maximum and minimum in each dimension of obtained Pareto set
if size(obtained_ps,1)==1
    obtained_min=obtained_ps;
    obtained_max=obtained_ps;
else
    obtained_min=min(obtained_ps);
    obtained_max=max(obtained_ps);
end
%find the maximum and minimum in each dimension of reference Pareto set
reference_min=min(reference_ps);
reference_max=max(reference_ps);
for i=1:n_var
    if reference_max(i)==reference_min(i)
        kesi(i)=1;% if the i_th value of reference Pareto set are all the same, ignore the i_th dimension.
    elseif obtained_min(i)>=reference_max(i)||reference_min(i)>=obtained_max(i)
         kesi(i)=0;
    else
        kesi(i)=((min(obtained_max(i),reference_max(i))-max(obtained_min(i),reference_min(i)))/...
        (reference_max(i)-reference_min(i)))^2;
    end
end
CR=nthroot(prod(kesi),2*n_var);
