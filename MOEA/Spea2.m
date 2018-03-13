function Spea2
  global SPEA2;

  SPEA2.defaultConfig = @defaultConfig_;
  SPEA2.run = @run_;
  SPEA2.name = "SPEA2";
end

function [result, h] = run_(population, ga_context, config)
  global SELECTION;
  global GA;
  
  N = config.N;
  M = config.M;
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
  
  archive = [];
  selection_fn = SELECTION.tournament(2);

  old_objective_values = [];
  g = 1;
  done = false;

  h_template = struct('population', [], 'objective_values', []);
  h = repmat(h_template, G_max, 1);
  
  while(~done)
	if (l == -1)
	  context.iteration = g;
    end
	
	%% Evaluation
	pool = vertcat(population, archive);
	[fitness, real_values_pop, objective_values] = evalFitnessAndPop_(pool, objective_vector, decode_fn, maximizing);
	indices_to_archive = environmentalSelection_(pool, objective_values, fitness, M);

	done = ((g == G_max) || stop_criteria_fn(objective_values, old_objective_values, maximizing));
	
	if (done)
	  non_dominated = fitness(indices_to_archive) < 1;
	  
	  %% TODO/FIXME: In case non_dominated is empty (no global non
	  %% dominated solution was found), we need to check if, among the
	  %% indices in archive, any individual is non dominated (we only
	  %% do a check on global elements (hence fitness < 1), but if no
	  %% global non dominated solution was found, some may be non
	  %% dominated _inside_ the solutions in the archive only).
	  %% For now, this asserts has not fired, but who knows...
      %%
	  %% assert(sum(non_dominated) > 0);
	  
      corresponding_values = real_values_pop(indices_to_archive, :);

	  result = corresponding_values(non_dominated, :);
	  
	  h(g).population = result;
	  h(g).objective_values = objective_values(non_dominated, :);
	else
	  archive = pool(indices_to_archive, :);
	  archive_fitness = fitness(indices_to_archive);

	  non_dominated = archive_fitness < 1;
	  h(g).population = archive(non_dominated, :);
	  h(g).objective_values = objective_values(non_dominated, :);
	  
	  %% Minus sign: tournament selection compare by >, but fitness is better when <.
	  selection = selection_fn(-archive_fitness);
	  mating_pool = archive(selection, :);

	  population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);
	  
	  old_objective_values = objective_values;
	  g = g + 1;
	end
  end
end

function result = defaultConfig_
  %DEFAULTCONFIG_ Preconfigured genetic algorithm config.
  %
  % Fields
  %  M                  Archive population count
  
  global GA;
  
  result = GA.defaultConfig();
  result.M = 20;
end

function [fitness, real_values_pop, objective_values] = evalFitnessAndPop_(pop_and_archive, fn_vector, decode_fn, maximizing)
  global UTILS;
  
  real_values_pop = decode_fn(pop_and_archive);
  objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);

  raw_fitness = rawFitness_(objective_values, maximizing);
  density_estimation = densityEstimation_(objective_values);
  fitness = raw_fitness + density_estimation;
end


function result = rawFitness_(objective_values, maximizing)
  [N, fn_count] = size(objective_values);
  BY_ROW = 2;

 result = zeros(1, N);
 strength = zeros(1, N);

 if (maximizing)
   relation = @ge; %% i dominates j, all objective values of i >= objective values of j
 else
   relation = @le;  %% i dominates j, all objective values of i <= objective values of j
 end
  
  %%TODO(@perf): Is there a way to remove the for-loop(s)?
  for i = 1:N
    %% dominated = all(objective_values(i, :) <= objective_values, BY_ROW); %% SLOWER!
	dominated = sum(relation(objective_values(i, :), objective_values), BY_ROW) == fn_count;
	strength(i) = sum(dominated) - 1; %% Because we also check the point itself.
  end
  
  for i = 1:N
    %% dominated = all(objective_values <= objective_values(i, :), BY_ROW); %% SLOWER!
	dominators = sum(relation(objective_values, objective_values(i, :)), BY_ROW) == fn_count;
	result(i) = sum(strength .* dominators') - strength(i); %% Because we also check the point itself.
  end
end

function result = densityEstimation_(objective_values)
  [N, ~] = size(objective_values);
  k = floor(sqrt(N));
  BY_ROW = 2;
  
  result = zeros(1, N);

  %%TODO(@perf): Is there a way to remove the for-loop?
  for i = 1:N
	distances = sum((objective_values - objective_values(i, :)).^2, BY_ROW);
	sorted_distances = mink(distances, k + 1); %% Because we also check the point itself.

	sigma = sorted_distances(end);  
	result(i) = 1 / (sigma + 2);
  end
end

function result = environmentalSelection_(pop_and_archive, objective_values, fitness, M)
  [N, ~] = size(pop_and_archive);
  
  non_dominated = (fitness < 1);
  non_dominated_count = sum(non_dominated);

  result = zeros(1, M);
  
  if (non_dominated_count == M)
	result  = non_dominated;
  elseif (non_dominated_count < M)
	indices = 1:N;
	result(1:non_dominated_count) = indices(non_dominated);

    fitness_rest = fitness;
	indices(non_dominated) = [];
    fitness_rest(non_dominated) = [];
    
    rest_to_fill = M - non_dominated_count;
    [~, sorted_indices] = mink(fitness_rest, rest_to_fill);
    
    can_fill_count = length(sorted_indices);
    
    %% NOTE: There may not even be enough individuals to fill the
    %% archive.
    %% So fill as much as it is possible.
    start_index = (non_dominated_count+1);
    end_index = (non_dominated_count+can_fill_count);
	result(start_index:end_index) = indices(sorted_indices);
    result(end_index+1:end) = [];
  else
	overcrowded_archive = pop_and_archive(non_dominated, :);
	overcrowded_archive_values = objective_values(non_dominated, :);
	
    valid_indices = find(non_dominated);
    indices_to_keep = truncationOperator_(overcrowded_archive, overcrowded_archive_values, M);
	result = valid_indices(indices_to_keep);
  end
end

function result = truncationOperator_(archive, objective_values, M) 
  [N, ~] = size(archive);
  BY_ROW = 2;

  to_remove_count = N - M;
  all_sorted_distances = zeros(N, N - 1);
  all_sorted_points = zeros(N, N - 1);
  
  for i = 1:N
	distances = sum((objective_values - objective_values(i, :)).^2, BY_ROW);
	[sorted_distances, sorted_points] = sort(distances);

	all_sorted_distances(i, :) = sorted_distances(2:end);  %% Because we also check the point itself.
	all_sorted_points(i, :) = sorted_points(2:end);
  end

  valid_points = 1:N;

  for i = 0:(to_remove_count-1)
    count = N - i;
    index_closest_point = comparePoints_(1:count, all_sorted_distances);
    point_to_remove = valid_points(index_closest_point);

    %% Remove point's info
    all_sorted_distances(index_closest_point, :) = [];
    all_sorted_points(index_closest_point, :) = [];
    
    indices_to_keep = (all_sorted_points ~= point_to_remove)';

	%% @perf
	%% @perf!
	%% @perf!!
	%% FIXME(@perf): This is the bottleneck of the program (and by far).
	%% This takes _way_ too much time (most is spent doing 'tempX(indices_to_keep)').
	%% And it grows in N^2!
    %% Find a way not to shift elements around...
    %% (Maybe a basic for loop... I do not know if it can be improved...)
	
    %% Remove, for all other points, their link to the removed point.
    %% (I need a temp because matlab does not allow "a'(foo)" syntax...)
    temp1 = all_sorted_distances';
    temp2 = all_sorted_points';
    
    all_sorted_distances = reshape(temp1(indices_to_keep), [], count - 1)';
    all_sorted_points = reshape(temp2(indices_to_keep), [], count - 1)';

    valid_points(index_closest_point) = [];
  end

  result = valid_points;
end

function result = comparePoints_(points, distances)
  [~, C] = size(distances);
  result = comparePointsRec_(points, distances, C);
end

function result = comparePointsRec_(points, distances, c)
  min_distance = min(distances(:, 1));
  corresponding_points = find(distances(:, 1) == min_distance);

  count_same = length(corresponding_points);
  
  if (count_same == 1)
	result = points(corresponding_points);
  elseif (c == 1)
	result = points(corresponding_points(randi(count_same, 1))); 
  else
    result = comparePointsRec_(corresponding_points, distances(corresponding_points, 2:end), c - 1);
  end
end
