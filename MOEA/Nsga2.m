function Nsga2
  global NSGA2;

  NSGA2.defaultConfig = @defaultConfig_;
  NSGA2.run = @run_;
end

function result = run_(population, ga_context, config)
  global SELECTION;
  global GA;
  
  N = config.N;
  l = config.l;
  
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
  
  selection_fn = SELECTION.tournament(2);

  old_objective_values = [];
  g = 1;
  
  [rank, ~, objective_values] = evalRankAndPop_(population, objective_vector, decode_fn, maximizing);

  %% First iteration, just select based on rank to create children.
  %% Minus sign: tournament selection compare by >, but rank is better when <.
  selection = selection_fn(-rank);
  parents = population;
  mating_pool = parents(selection, :);
  
  population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);

  old_objective_values = objective_values;
  g = g + 1;

  while (~((g == G_max) || stop_criteria_fn(objective_values, old_objective_values, maximizing)))
	pool = vertcat(parents, population);
	
	%% Evaluation
	[rank, ~, objective_values] = evalRankAndPop_(pool, objective_vector, decode_fn, maximizing);

	%% TODO: Separate into own function?
	new_population_count = 0;
	indices_to_keep = zeros(1, N);

	fitness_indices_to_keep = zeros(1, N);

	no_overfill = true;
	front_index = 1;

	%% There is no do-while...
	while (no_overfill)
	  belong_to_front = (rank == front_index);
	  front_indices = find(belong_to_front);
	  
	  belong_to_front_count = length(front_indices);

	  front_crowding_distance = crowdingDistanceAssignment_(objective_values(front_indices, :));

	  %% As long as we do not overfill, we can add all of it
	  if ((belong_to_front_count + new_population_count) <= N)
		append_indices = (new_population_count + 1):(new_population_count + belong_to_front_count);
		
		indices_to_keep(append_indices) = front_indices;
		
		[~, sorted_indices] = sort(-front_crowding_distance);
        bias = zeros(1, belong_to_front_count);
		bias(sorted_indices) = (1:belong_to_front_count) / (belong_to_front_count + 1);
		fitness_indices_to_keep(append_indices) = front_index + bias;

		new_population_count = new_population_count + belong_to_front_count;
        
        if (new_population_count == N)
          break
        else
          front_index = front_index + 1;
        end
	  else
		no_overfill = false;
	  end
	end

	%% We can partially add the last front
	if (new_population_count < N)
	  remaining_count = N - new_population_count;
	  
	  [~, sorted_indices] = mink(-front_crowding_distance, remaining_count);

	  abs_indices = front_indices(sorted_indices);
	  
	  indices_to_fill = (new_population_count+1):N;
	  indices_to_keep(indices_to_fill) = abs_indices;
	  
	  bias = (1:remaining_count) / (remaining_count + 1);
	  fitness_indices_to_keep(indices_to_fill) = front_index + bias;
    end

	parents = pool(indices_to_keep, :);
    selection = selection_fn(-fitness_indices_to_keep);
    mating_pool = parents(selection, :);

	population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);

	old_objective_values = objective_values;
	g = g + 1;
  end

  [rank, real_values_pop, ~] = evalRankAndPop_(population, objective_vector, decode_fn, maximizing);
  
  result = real_values_pop((rank == 1), :);
end

function result = defaultConfig_
	 %DEFAULTCONFIG_ Preconfigured genetic algorithm config.
	 %
	 % Fields
	 %  N                  Population count
	 %  G_max              Max iteration count
	 %  l                  Chromosome length, in [1, 53]
	 %  Pc                 Crossover probability
	 %  Pm                 Mutation probability
	 %  crossover_fn       Crossover function
	 %  mutation_fn        Mutation function
	 %  stop_criteria_fn   Stop criteria function
	 %  clamp_fn           Clamp function, not used with binary values
	 %
	 % See also Selection, Crossover, Mutation,
	 % StopCriteria, Clamp.
  
  global CROSSOVER;
  global MUTATION;
  global STOP_CRITERIA;
  global CLAMP;
  
  result.N = 100;
  result.G_max = 100;
  
  %% NOTE: 'binary' is just an integer representation (to get to the
  % actual value => v = (i / maxI) * (c(1) - c(0)) + c(0), with c the
  % constaints for this variable)
  result.l = 12;
  
  result.Pc = 0.5;
  result.Pm = 0.1;

  result.crossover_fn = CROSSOVER.singlePoint;
  result.mutation_fn = MUTATION.bitFlip;
  result.stop_criteria_fn = STOP_CRITERIA.time;
  result.clamp_fn = CLAMP.default;
end

function [rank, real_values_pop, objective_values] = evalRankAndPop_(population, fn_vector, decode_fn, maximizing)
  global UTILS;
  
  real_values_pop = decode_fn(population);
  objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);

  rank = fastNonDominatedSort_(objective_values, maximizing);
end

function result = fastNonDominatedSort_(objective_values, maximizing)
  [N, fn_count] = size(objective_values);
  BY_COLUMN = 1;
  BY_ROW = 2;

  result = zeros(1, N);
  
  if (maximizing)
    relation = @ge; %% i dominates j, all objective values of i >= objective values of j
  else
    relation = @le;  %% i dominates j, all objective values of i <= objective values of j
  end

  point_domination = zeros(N, N);
  
  for i = 1:N
      domination = sum(relation(objective_values(i, :), objective_values), BY_ROW)' == fn_count;
    
      %% In case multiple points have the same values, we are both dominated and dominators of those other points.
      %% In that case, we do not want to consider them as being dominated or dominators.
      %% As we also check the point itself, this also takes care of that.
      is_dominated = (domination == 1);
      is_same = sum(objective_values(i, :) == objective_values(is_dominated, :), BY_ROW)' == fn_count;
      
      domination(is_dominated) = ~is_same;
      point_domination(i, :) = domination;
  end
  
  point_dominators_count = sum(point_domination, BY_COLUMN);

  indices_to_compare = 1:N;
  indices_to_compare_count = N;
  
  front_index = 1;
  while (indices_to_compare_count)
	rel_front_indices = find(point_dominators_count(indices_to_compare) == 0);
	front_indices_count = length(rel_front_indices);
	
	assert(front_indices_count ~= 0);

	abs_front_indices = indices_to_compare(rel_front_indices);
	result(abs_front_indices) = front_index;

	point_dominators_count = point_dominators_count - sum(point_domination(abs_front_indices, :), BY_COLUMN);

	front_index = front_index + 1;

	indices_to_compare(rel_front_indices) = [];
	indices_to_compare_count = indices_to_compare_count - front_indices_count;
  end
end

function result = crowdingDistanceAssignment_(objective_values)
  [N, fn_count] = size(objective_values);
  result = zeros(1, N);

  for m = 1:fn_count
	values = objective_values(:, m);
	
	[~, indices] = sort(values);
	first_index = indices(1);
	last_index = indices(end);

	%% Max - min
	values_interval = values(last_index) - values(first_index);
	
	result(first_index) = Inf;
	result(last_index) = Inf;

	for i = 2:(N-1)
	  result(i) = result(i) + (values(i + 1) - values(i - 1)) / values_interval;
	end
  end
end
