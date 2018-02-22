function Spea2
  global SPEA2;

  SPEA2.optimize = @optimize;
  SPEA2.maximize = @maximize;
  SPEA2.minimize = @minimize;
  
  SPEA2.showPaleto = @showPaleto_;

  SPEA2.defaultConfig = @defaultConfig;
end

function result = defaultConfig
	 %DEFAULTCONFIG Preconfigured genetic algorithm config.
	 %
	 % Fields
	 %  N                  Population count
	 %  M                  Archive population count
	 %  G_max              Max iteration count
	 %  l                  Chromosome length, in [1, 53]
	 %  Pc                 Crossover probability
	 %  Pm                 Mutation probability
	 %  selection_fn       Selection function
	 %  crossover_fn       Crossover function
	 %  mutation_fn        Mutation function
	 %  stop_criteria_fn   Stop criteria function
	 %  clamp_fn           Clamp function, not used with binary values
	 %
	 % See also Selection, Crossover, Mutation,
	 % StopCriteria, Clamp.
  
  global SELECTION;
  global CROSSOVER;
  global MUTATION;
  global STOP_CRITERIA;
  global CLAMP;
  
  result.N = 100;
  result.M = 20;
  result.G_max = 100;
  
  %% NOTE: 'binary' is just an integer representation (to get to the
  % actual value => v = (i / maxI) * (c(1) - c(0)) + c(0), with c the
  % constaints for this variable)
  result.l = 12;
  
  result.Pc = 0.5;
  result.Pm = 0.1;

  result.selection_fn = SELECTION.wheel;
  result.crossover_fn = CROSSOVER.singlePoint;
  result.mutation_fn = MUTATION.bitFlip;
  result.stop_criteria_fn = STOP_CRITERIA.time;
  result.clamp_fn = CLAMP.default;
end

function result = optimize(maximizing, objective_vector, constraints, config)
  global UTILS;
  global SELECTION;

  %% TODO: Parameter check and default value.
  %% TODO: Make sure only binary crossover functions can be used if
  %% l >= 1.
  %% Same for arithmetic functions and l == -1.
  
  N = config.N;
  M = config.M;
  l = config.l;
  
  G_max = config.G_max;

  Pc = config.Pc;
  Pm = config.Pm;

  selection_fn = SELECTION.tournament(2, N);
  crossover_fn = config.crossover_fn;
  mutation_fn = config.mutation_fn;
  stop_criteria_fn = config.stop_criteria_fn;
  clamp_fn = config.clamp_fn;

  decode_fn = UTILS.decode(constraints, l);
  
  if (l == -1)
	context = struct('constraints', constraints, 'G_max', G_max, 'iteration', 0, 'clamp_fn', clamp_fn);
  else
	context = l;
  end
  
  tic;
  
  [var_count, ~] = size(constraints);

  population = initialGeneration_(N, constraints, l);
  archive = [];
  
  last_iteration = G_max + 1;
  old_fitness = [];
  g = 1;
  while (g <= G_max)
	if (l == -1)
	  context.iteration = g;
	end

	%% Evaluation
	%% TODO: Should we add a fitness_vector parameter and use it here
	%% instead? (Like we did with Ga)
	pool = vertcat(population, archive);
	[fitness, real_values_pop, objective_values] = evalFitnessAndPop_(pool, objective_vector, decode_fn, maximizing);
	indices_to_archive = environmentalSelection_(pool, objective_values, fitness, M);

	%%TODO: Stop criteria on objective_values / fitness_values (if we
	%% change objective_vector to fitness_vector).
	%% Change StopCriteria to take a cell array of functions, and a
	%% matrix of results. (not just one function and an array)
	%% if ((g == G_max) || stop_criteria_fn(fitness, old_fitness))
	if (g == G_max)
	  last_iteration = g;

	  non_dominated = fitness(indices_to_archive) < 1;
	  
	  %% TODO/FIXME: In case non_dominated is empty (no global non
	  %% dominated solution was found), we need to check if, among the
	  %% indices in archive, any individual is non dominated (we only
	  %% do a check on global elements (hence fitness < 1), but if no
	  %% global non dominated solution was found, some may be non
	  %% dominated _inside_ the solutions in the archive only).
	  %% For now, this asserts has not fired, but who knows...
	  assert(~isempty(non_dominated));
	  
      corresponding_values = real_values_pop(indices_to_archive, :);

	  result = corresponding_values(non_dominated, :);
	  break;
    end
	
	archive = pool(indices_to_archive, :);
    archive_fitness = fitness(indices_to_archive);

    %% Minus sign: tournament selection compare by >, but fitness is better when <.
    %% TODO: Make tournament selection work with <, or keep negation...
    selection = selection_fn(-archive_fitness);
    mating_pool = archive(selection, :);

	children = crossover_(mating_pool, crossover_fn, Pc, context);

	%% Every allele that needs to mutate is 1 at the correponding index
	if (l == -1)
	  mutations = rand(N, var_count, 1) <= Pm;
	else
	  mutations = rand(N, l, var_count) <= Pm;
	end

	population = mutation_fn(children, mutations, context);


	old_fitness = fitness;
	g = g + 1;
  end

  toc;
end

function result = maximize(objective_vector, constraints, config)
  result = optimize(1, objective_vector, constraints, config);
end

function result = minimize(objective_vector, constraints, config)
  result = optimize(0, objective_vector, constraints, config);
end

function [fitness, real_values_pop, objective_values] = evalFitnessAndPop_(pop_and_archive, fn_vector, decode_fn, maximizing)
  global UTILS;
  
  real_values_pop = decode_fn(pop_and_archive);
  objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);

  raw_fitness = rawFitness_(objective_values, maximizing);
  density_estimation = densityEstimation_(objective_values);
  fitness = raw_fitness + density_estimation;
end

function result = initialGeneration_(N, constraints, l)
  global UTILS;
  
  if (l == -1)
	result = UTILS.randomIn(constraints, N);
  else
	dim = size(constraints);
	var_count = dim(1);

	max_val = 2^l-1;

	result = randi(max_val, N, var_count);
  end
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
  global UTILS;
  
  [N, var_count] = size(pop_and_archive);
  
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
    
    %% NOTE: There may not be even enough individuals to fill the
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
  %[~, C] = size(distances);
  min_distance = min(distances(:, 1));
  corresponding_points = find(distances == min_distance);

  count_same = length(corresponding_points);
  
  if (count_same == 1)
	result = points(corresponding_points);
  %elseif (C == 1)
  %  result = points(corresponding_points(randi(count_same, 1))); 
  else
    result = comparePoints_(corresponding_points, distances(corresponding_points, 2:end));
  end
end

function result = crossover_(mating_pool, crossover_fn, Pc, context)
  
  %% Modify mating pool to have an array of [i, j] (two individuals on
  %% the same row), so we do not have to introduce an explicit loop
  %% (usually slower) to compute the crossover of each parent pair.
  var_count = length(mating_pool(1, :));
  mating_pool = reshape(mating_pool', 2 * var_count, [])';

  rand_val = rand(length(mating_pool(:, 1)), 1);
  indices = find(rand_val <= Pc); %% Find which pair will crossover

  go_through_crossover = mating_pool(indices, :);

  %% Pair separation
  min_b = 1:var_count; %% The first half (first individual)
  max_b = (var_count+1):(var_count * 2); %% The second half (second individual)

  %% NOTE/TODO?: Instead of separating the variables, we could
  %% concatenated them (shift each by l*i bits and directly apply the
  %% crossover and mutation over the resulting integer) and split them
  %% after everything is done. This could, potentialy, speed up the
  %% computation.
  %% Howewer, octave has a limit of 53 bits (at least, bitset limits
  %% the index to 53), which means we would be limited to 53 /
  %% var_count bits per variable. (var_count = 3 => 17 bits)
  %% (By handling them separately, we do not have _any_ limitation)
  unchanged = mating_pool;
  unchanged(indices, :) = [];  %% Remove pairs which are going to crossover.

  go_through_crossover = crossover_fn(go_through_crossover(:, min_b), go_through_crossover(:, max_b), context);

  %% Flatten the result to have [i1; i2; ...] again, instead of
  %% [ [i1, i2]; [i3, i4]; ... ]
  result = reshape([unchanged; go_through_crossover]', var_count, [])';
end

%% TODO: Test 3D function plot
%%       Add legend
function showPaleto_(problem, variables)
  global UTILS;
  figure(5)
  clf;

  constraints = problem.constraints;
  objective_vector = problem.objective_vector;

  [~, var_count] = size(variables);
  [~, fn_count] = size(objective_vector);

  can_plot_var = (var_count <= 3);
  can_plot_fn = (fn_count <= 3);
  plot_count = can_plot_var + can_plot_fn;
  
  if (var_count == 2)
      m = 50;
  else
      m = 10;
  end
  
  if (var_count == 2)
    [dx, dy] = meshgrid(UTILS.linspacea(constraints, m));
    dz = 0 * dx;
  else
    [dx, dy, dz] = meshgrid(UTILS.linspacea(constraints, m));
  end
  
  if (~can_plot_var)
	warning('showPaleto: variables: can not plot in 4D, sorry!');
  else
	subplot(1, plot_count, 1);
    hold on;
    
    %% TODO: Plot domain in 3D (plot3(dx, dy, dz) does not work)
    if (var_count == 3)
        %%mesh(dx, dy, dz);
        plotN_(variables, @plot3, 'r+');  %% TODO: Test
    else    
        mesh(dx, dy, dz);
        plotN_(variables, @plot, 'r+');
    end
  end

  if (~can_plot_fn)
	warning('showPaleto: variables: can not plot in 4D, sorry!');
  else
	subplot(1, plot_count, plot_count);
    hold on;

    if (var_count == 2)
        fz = evalFn_(objective_vector, dx, dy);
    else
        fz = evalFn_(objective_vector, dx, dy, dz);
    end
    
	values = UTILS.evalFnVector(objective_vector, variables);
    
    %% TODO: Plot domain in 3D (no 3D function yet)
	if (fn_count == 3)
        mesh(fz(:, :, 1), fz(:, :, 2), fz(:, :, 3));
        plotN_(values, @plot3, 'r+'); %% TODO: Test
    else
        x = fz(:, :, 1);
        y = fz(:, :, 2);
        z = zeros(size(x));
        
        C = gradient(x .* y);
        mesh(x, y, z, C);
        plotN_(values, @plot, 'r+');
    end
  end
end

function plotN_(val_array, plot_fn, style)
    BY_COLUMN = 2; 
    to_var_arg = num2cell(val_array', BY_COLUMN);
    
    plot_fn(to_var_arg{:}, style)
end

%% NOTE: This is needed to allow evaluation of objective_vector on
%% matrices (which is not handle with UTILS.evalFn (which was not
%% designed to handle it anyway).
%% Btw, this is an ugly hack!)
function result = evalFn_(objective_vector, varargin)
  [~, fn_count] = size(objective_vector);

  %% It works. I do reall know why, by it does!
  [N, M] = size(varargin{1});
  result = zeros(M, N, fn_count);

  for i = 1:fn_count
	f = objective_vector{i};

	result(:, :, i) = f(varargin{:});
  end
end
