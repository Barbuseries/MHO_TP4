function Pesa2
  global PESA2;

  PESA2 = Pesa_common("Pesa2", @region_based_selection_);
end

function h = region_based_selection_(population_grid_index, all_grid_indices, all_squeeze_factors, selection_fn)
  h = @(count) get_individual_by_region_(population_grid_index, all_grid_indices(selection_fn(-all_squeeze_factors, count)));
end

function result = get_individual_by_region_(population_grid_index, grid_indices)
  c = length(grid_indices);
  
  result = zeros(1, c);
  
  %% TODO (@perf): This can be improved by, for a given grid index,
  %% only querying individuals with the same index once (and count the
  %% number of times grid index appears, so we just select K
  %% individuals once instead of selecting one individual K times).
  for i = 1:c
	possible_individuals = find(population_grid_index == grid_indices(i));
    
	result(i) = possible_individuals(randi(length(possible_individuals), 1));
  end
end
