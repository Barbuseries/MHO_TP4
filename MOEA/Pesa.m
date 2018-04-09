function Pesa
  global PESA;

  PESA.defaultConfig = @defaultConfig_;
  PESA.run = @run_;
  PESA.name = "PESA";
end

function [result, h] = run_(internal_population, ga_context, config)
  global SELECTION;
  global GA;
  
  N = config.N; %% Internal population size
  M = config.M; %% External population size
  l = config.l;
  C = config.C;
  
  G_max = config.G_max;

  Pc = config.Pc;
  Pm = config.Pm;

  crossover_fn = config.crossover_fn;
  mutation_fn = config.mutation_fn;
  stop_criteria_fn = config.stop_criteria_fn;
    
  objective_vector = ga_context.objective_vector;
  maximizing = ga_context.maximizing;
  decode_fn = ga_context.decode_fn;
  context = ga_context.operator_context;

  external_population = [];
  external_population_objective_values = [];
  external_population_grid_index = [];
  
  grid_indices = [];
  grid_squeeze_factors = [];
  
  selection_fn = SELECTION.tournament(2);

  hyper_grid_fn = hyper_grid_indexing_(ga_context.constraints, C);

  old_external_population_objective_values = [];
  g = 1;
  done = false;

  h_template = struct('population', [], 'objective_values', []);
  h = repmat(h_template, G_max, 1);
  
  while(~done)
	if (l == -1)
	  context.iteration = g;
    end
	
	%% Evaluation
	[non_dominated, real_values_pop, objective_values] = evalFitnessAndPop_(internal_population, objective_vector, decode_fn, maximizing);

	population_to_archive = internal_population(non_dominated, :);
	corresponding_objective_values = objective_values(non_dominated, :);
	corresponding_grid_index = hyper_grid_fn(real_values_pop(non_dominated, :));
	
	[external_population, external_population_objective_values, external_population_grid_index, squeeze_factor, grid_indices, grid_squeeze_factors] = ...
	updateArchive_(population_to_archive, corresponding_objective_values, corresponding_grid_index, ...
				   external_population, external_population_objective_values, external_population_grid_index, maximizing, M, grid_indices, grid_squeeze_factors);

	done = ((g == G_max) || stop_criteria_fn(external_population_objective_values, old_external_population_objective_values, maximizing));
	
	if (done)
	  result = decode_fn(external_population);
	  
	  h(g).population = result;
	  h(g).objective_values = external_population_objective_values;
	else
	  crossover_count = sum(rand(N, 1) <= Pc); %% Crossovers happen with probability Pc
	  mutation_count  = N - crossover_count; %% Otherwhise, it's a mutation

	  crossover_selection = selection_fn(-squeeze_factor, 2 * crossover_count);
	  mutation_selection  = selection_fn(-squeeze_factor, mutation_count);

	  crossover_pool = external_population(crossover_selection, :);
	  mutation_pool  = external_population(mutation_selection, :);

	  children_from_crossover = GA.crossover(crossover_pool, crossover_fn, context);
	  children_from_mutation = GA.mutate(mutation_pool, l, mutation_fn, Pm, context);

	  %% Keep one child out of the two children generated (with
	  %% probablity = 1/2).
	  child_to_keep = rand(crossover_count, 1) <= 0.5;
	  first_child = (0:2:((2*crossover_count)-1)) .* child_to_keep';
	  second_child = (1:2:((2*crossover_count)-1)) .* ~child_to_keep';

      child_to_keep_indices = (first_child + second_child + 1)';
	  children_kept_from_crossover = children_from_crossover(child_to_keep_indices, :);
	  
	  internal_population = [children_kept_from_crossover; children_from_mutation];


	  h(g).population = decode_fn(external_population);
	  h(g).objective_values = external_population_objective_values;
	  
	  old_external_population_objective_values = external_population_objective_values;
	  g = g + 1;
	end
  end
end

function result = defaultConfig_
  %DEFAULTCONFIG_ Preconfigured genetic algorithm config.
  %
  % Fields
  %  M                  External population count
  %  C                  Size of the hyper-grid (CxC)
  
  global GA;
  
  result = GA.defaultConfig();
  result.M = 20;
  result.C = 32;
end

function [non_dominated, real_values_pop, objective_values] = evalFitnessAndPop_(population, fn_vector, decode_fn, maximizing)
  global UTILS;
  
  real_values_pop = decode_fn(population);
  objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);

  is_dominated = domination_(objective_values, maximizing);
  non_dominated = find(is_dominated == 0);
end

function result = domination_(objective_values, maximizing)
  [N, fn_count] = size(objective_values);
  BY_ROW = 2;

  result = zeros(1, N);

  if (maximizing)
	relation = @ge; %% i dominates j, all objective values of i >= objective values of j
  else
	relation = @le;  %% i dominates j, all objective values of i <= objective values of j
  end
  
  for i = 1:N
	dominators = sum(relation(objective_values, objective_values(i, :)), BY_ROW) == fn_count;
    is_equal = sum(objective_values == objective_values(i, :), BY_ROW) == fn_count;
    
	result(i) = sum(dominators & ~is_equal) ~= 0;
  end
end

function [dominators, dominates] = domination_status_(a, b, maximizing)  
  [N, fn_count] = size(a);
  [M, ~] = size(b);
  
  BY_ROW = 2;

  dominates = zeros(N, M);
  dominators = zeros(N, M);
  
  if (isempty(b))
      return;
  end

  if (maximizing)
	relation = @ge; %% i dominates j, all objective values of i >= objective values of j
  else
	relation = @le;  %% i dominates j, all objective values of i <= objective values of j
  end

  domination_fn = @(x, y) (sum(relation(x, y), BY_ROW) == fn_count) & ~(sum(x == y, BY_ROW) == fn_count);
  
  for i = 1:N
	dominates(i, :) = domination_fn(a(i, :), b);
  end
  
  for i = 1:N
	dominators(i, :) = domination_fn(b, a(i, :));
  end
end

function [external_population, external_population_objective_values, ...
          external_population_grid_index, squeeze_factor, ...
          all_grid_indices, all_squeeze_factors] = updateArchive_(population, population_objective_values, ...
                                                                  population_grid_index, external_population, ...
                                                                  external_population_objective_values, ...
                                                                  external_population_grid_index, maximizing, M, ...
                                                                  all_grid_indices, all_squeeze_factors)
  [population_count, ~] = size(population);
  [old_M, ~] = size(external_population);
  
  current_M = old_M;

  %% NOTE: Because we know population only contains non-dominated
  %% individuals, there is no need to, 1: check the domination within
  %% the population, 2: update the domination status when an
  %% individual is added to the archive (because we know it will be
  %% added _at the end_ and that we do not remove any element in the
  %% archive during the loop below).
  %% (This is done by both setting the grid index of a removed
  %% individual to -1, and their column in the domination status to 0
  %% (NOTE: Therefore only archive members present at the start of the
  %% loop are indexed in this table)).
  [dominators, dominates] = domination_status_(population_objective_values, external_population_objective_values, maximizing);

  for i = 1:population_count
	individual_dominators = dominators(i, :);

	if (sum(individual_dominators) ~= 0)
	  continue;
	end
	
	grid_index = population_grid_index(i);

	archive_dominated = find(dominates(i, :));
	archive_dominated_count = length(archive_dominated);

	%% Get, for a given grid index, the number of individuals to remove.
	[grid_index_to_update, squeeze_factor_diff] = get_squeeze_factor_by_index_(external_population_grid_index(archive_dominated));

	%% Individual is now inactive and will be removed after the loop is done.
	external_population_grid_index(archive_dominated) = -1;

	%% An inactive element is not relevent anymore as a dominator or
	%% as dominated.
	dominators(:, archive_dominated) = 0;
	dominates(:, archive_dominated) = 0;

	modified_squeeze_factor_count = length(grid_index_to_update);

	updated_new_element_grid_index_during_remove = 0;
	start_index_search = 1;
	for j = 1:modified_squeeze_factor_count
	  update_index = grid_index_to_update(j);
	  update_diff = squeeze_factor_diff(j);

	  %% NOTE: As all_grid_indices is guaranteed to be sorted, we can
	  %% start the search after the previously found index.
	  index_in_table = (start_index_search - 1) + find(all_grid_indices(start_index_search:end) == update_index, 1);
	  entry_squeeze_factor = all_squeeze_factors(index_in_table);

	  index_is_same_as_new_element = (update_index == grid_index);
	  all_squeeze_factors(index_in_table) = entry_squeeze_factor - update_diff + index_is_same_as_new_element;

	  if (index_is_same_as_new_element)
		updated_new_element_grid_index_during_remove = 1;
	  end

	  start_index_search = index_in_table + 1;
	end

	if (~updated_new_element_grid_index_during_remove)
	  index_in_table = find(all_grid_indices >= grid_index, 1);

	  if (isempty(index_in_table))
		all_grid_indices = [all_grid_indices, grid_index];
		all_squeeze_factors = [all_squeeze_factors, 1];
      elseif (all_grid_indices(index_in_table) > grid_index)
		all_grid_indices(index_in_table:end+1) = [grid_index, all_grid_indices(index_in_table:end)];
		all_squeeze_factors(index_in_table:end+1) = [1, all_squeeze_factors(index_in_table:end)];
      else
          all_squeeze_factors(index_in_table) = all_squeeze_factors(index_in_table) + 1;
	  end
	end

	%% TODO (@perf): See if this takes too much time.
	%% (Instead, we could pre-resize the external population to be N +
	%% M and just strip the excess later).
	external_population = [external_population; population(i, :)];
	external_population_objective_values = [external_population_objective_values; population_objective_values(i, :)];
	external_population_grid_index = [external_population_grid_index, population_grid_index(i)];

	current_M = current_M - archive_dominated_count + 1;

	if (current_M > M)
	  [max_squeeze_factor, ~] = max(all_squeeze_factors);
	  all_max_squeeze_factors_indices = find(all_squeeze_factors == max_squeeze_factor);

      squeeze_factor_to_update = all_max_squeeze_factors_indices(randi(length(all_max_squeeze_factors_indices), 1));
	  grid_index_to_update = all_grid_indices(squeeze_factor_to_update);
	  associated_indices_in_pop = find(external_population_grid_index == grid_index_to_update);
	  
	  random_individual_to_remove = associated_indices_in_pop(randi(length(associated_indices_in_pop), 1));
	  
	  all_squeeze_factors(squeeze_factor_to_update) = all_squeeze_factors(squeeze_factor_to_update) - 1;

	  %% The individual was present in the original external
	  %% population, therefore we need to update the domination
	  %% status (to not take it into account anymore).
	  if (random_individual_to_remove <= old_M)
		dominators(:, random_individual_to_remove) = 0;
		dominates(:, random_individual_to_remove) = 0;
	  end
	  
	  external_population_grid_index(random_individual_to_remove) = -1;
      current_M = current_M - 1;
	end
  end

  grid_indices_to_remove = find(all_squeeze_factors == 0);
  all_grid_indices(grid_indices_to_remove) = [];
  all_squeeze_factors(grid_indices_to_remove) = [];

  %% Remove inactive archive elements
  external_population_to_remove = find(external_population_grid_index == -1);
  external_population(external_population_to_remove, :) = [];
  external_population_objective_values(external_population_to_remove, :) = [];
  external_population_grid_index(external_population_to_remove) = [];

  %% Update squeeze factors
  [i, ~] = find(external_population_grid_index == all_grid_indices');
  squeeze_factor = all_squeeze_factors(i);
end

function h = hyper_grid_indexing_(constraints, C)
  [var_count, ~] = size(constraints);
  
  min_values = constraints(:, 1)';
  max_values = constraints(:, 2)';
  
  %% +1: Avoid having h(max_value) == (C + 1).
  grid_width = (max_values + 1 - min_values) / C;

  powers = C.^(0:(var_count - 1));

  h = @(real_pop) hyper_grid_indexing_inner_(real_pop, min_values, grid_width, powers);
end

function result = hyper_grid_indexing_inner_(real_values_pop, min_values, grid_width, powers)
  BY_ROW = 2;
  
  result = sum(floor((real_values_pop - min_values) / grid_width) .* powers, BY_ROW)';
end

function [all_indices, all_squeeze_factors] = get_squeeze_factor_by_index_(grid_index)
  if (isempty(grid_index))
      all_indices = [];
      all_squeeze_factors = [];
  else
      [a, ~, c] = unique(grid_index);
      result = [a', accumarray(c, 1)];

      all_indices = result(:, 1)';
      all_squeeze_factors = result(:, 2)';
  end
end
