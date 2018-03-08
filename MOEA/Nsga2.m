function Nsga2
  global NSGA2;

  NSGA2.defaultConfig = @defaultConfig_;
  NSGA2.init = @init_;
  NSGA2.step = @step_;
end

function result = init_(config)
  global SELECTION;
  
  result = struct('parents', [], ...
				  'selection_fn', SELECTION.tournament(2, config.N));
end

function [done, result, objective_values, my_context] = step_(population, ga_context, my_context)
  [N, ~] = size(population);
  objective_vector = ga_context.objective_vector;
  maximizing = ga_context.maximizing;
  
  stop_criteria_fn = ga_context.stop_criteria_fn;
  decode_fn = ga_context.decode_fn;

  g = ga_context.iteration;
  G_max = ga_context.G_max;
  old_objective_values = ga_context.old_objective_values;
  
  selection_fn = my_context.selection_fn;
  parents = my_context.parents;
  
  done = false;

  pool = vertcat(parents, population);
  
  %% Evaluation
  [rank, ~, objective_values] = evalRankAndPop_(pool, objective_vector, decode_fn, maximizing);

  %% First iteration, just select based on rank to create children.
  if (g == 1)
	%% Minus sign: tournament selection compare by >, but fitness is better when <.
	selection = selection_fn(-rank);
	result = population(selection, :);
  else
	new_population_count = 0;
	indices_to_keep = zeros(1, N);

	no_overfill = true;
	front_index = 1;

	%% There is no do-while...
	while (no_overfill)
	  belong_to_front = (rank == front_index);
	  belong_to_front_count = sum(belong_to_front);

	  %% As long as we do not overfill, we can add all of it
	  if ((belong_to_front_count + new_population_count) <= N)
		append_indices = (new_population_count + 1):(new_population_count + belong_to_front_count);
		indices_to_keep(append_indices) = find(belong_to_front);

		new_population_count = new_population_count + belong_to_front_count;
        
        front_index = front_index + 1;
	  else
		no_overfill = false;
      end
	end

	%% We can partially add the last front
	if (population_count < N)
	  %%TODO: Crowding distance sort on front.
	  error('Not yet done, sorry');
    end
    
    result = pool(indices_to_keep, :);
  end

  if ((g == G_max) || stop_criteria_fn(objective_values, old_objective_values, maximizing))
	done = true;
    
    result = decode_fn(result);
  end

  my_context.parents = population;
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

%% NOTE: Some liberties where taken to implement this function. It
%% should return the same result as the original one, but this one
%% should be faster (tests incoming!).
%% This was done because this version is easier to write and to
%% vectorise (so Matlab can run it quicker).
function result = fastNonDominatedSort_(objective_values, maximizing)
  [N, fn_count] = size(objective_values);
  BY_ROW = 2;

  result = zeros(1, N);
  indices_to_compare = 1:N;
  indices_to_compare_count = N;
  
  if (maximizing)
    relation = @gt; %% i is dominated by j, all objective values of i > objective values of j
  else
    relation = @lt;  %% i is dominated by j, all objective values of i < objective values of j
  end

  front_index = 1;
  while (indices_to_compare_count > 0)
	points_to_compare = objective_values(indices_to_compare, :);
	is_non_dominated = zeros(1, indices_to_compare_count);

	for i = 1:indices_to_compare_count
	  is_dominated = (sum(relation(points_to_compare, points_to_compare(i, :)), BY_ROW) == fn_count);
	  is_non_dominated(i) = (sum(is_dominated) == 0);
    end
    
	rel_non_dominated_indices = is_non_dominated ~= 0;
    abs_non_dominated_indices = indices_to_compare(rel_non_dominated_indices);
    
	result(abs_non_dominated_indices) = front_index;

    indices_to_compare(rel_non_dominated_indices) = [];
	indices_to_compare_count = indices_to_compare_count - length(abs_non_dominated_indices);
    
	front_index = front_index + 1;
  end
end
