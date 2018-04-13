function Ga
  global GA;

  GA.iterativeShowPaleto = @iterativeShowPaleto_;
  GA.showPaleto = @showPaleto_;
  GA.defaultConfig = @defaultConfig_;
  GA.make_new_pop = @make_new_pop;
  GA.initialPopulation = @initialPopulation_;
  GA.crossover = @crossover_all_;
  GA.mutate = @mutate_;
  
  GA.create = @create_ga_;
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

function [result, h] = optimize(ga, maximizing, objective_vector, constraints, config)
  global UTILS;
  
  %% TODO: Parameter check and default value.
  %% TODO: Make sure only binary crossover functions can be used if
  %% l >= 1.
  %% Same for arithmetic functions and l == -1.

  N = config.N;
  l = config.l;
  
  G_max = config.G_max;

  decode_fn = UTILS.decode(constraints, l);

  if (l == -1)
	context = struct('constraints', constraints, 'G_max', G_max, 'iteration', 0, 'clamp_fn', config.clamp_fn);
  else
	context = l;
  end

  ga_context = struct('objective_vector', {objective_vector}, ...
					  'constraints', constraints, ...
					  'maximizing', maximizing, ...
					  'decode_fn', decode_fn, ...
					  'operator_context', context);
  
  %%tic;

  if (isfield(config, 'population') && ~isempty(config.population))
	population = config.population;
  else
	population = initialPopulation_(N, constraints, l);
  end
  
  [result, h] = ga.run(population, ga_context, config);
  
  %%toc;
end

function [result, h] = maximize(ga, objective_vector, constraints, config)
  [result, h] = optimize(ga, 1, objective_vector, constraints, config);
end

function [result, h] = minimize(ga, objective_vector, constraints, config)
  [result, h] = optimize(ga, 0, objective_vector, constraints, config);
end

function result = initialPopulation_(N, constraints, l)
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

function result = create_ga_(ga)
  result = ga;
  
  result.optimize = @(varargin) optimize(ga, varargin{:});
  result.maximize = @(varargin) maximize(ga, varargin{:});
  result.minimize = @(varargin) minimize(ga, varargin{:});
end

function result = crossover_all_(mating_pool, crossover_fn, context)
  [N, var_count] = size(mating_pool);
  
  if (N == 0)
      result = [];
      return;
  end

  %% Modify mating pool to have an array of [i, j] (two individuals on
  %% the same row), so we do not have to introduce an explicit loop
  %% (usually slower) to compute the crossover of each parent pair.
  mating_pool = reshape(mating_pool', 2 * var_count, [])';

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
  after_crossover = crossover_fn(mating_pool(:, min_b), mating_pool(:, max_b), context);

  %% Flatten the result to have [i1; i2; ...] again, instead of
  %% [ [i1, i2]; [i3, i4]; ... ]
  result = reshape(after_crossover', var_count, [])';
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

function result = mutate_(population, l, mutation_fn, Pm, context)
  [N, var_count] = size(population);
  
  %% Every allele that needs to mutate is 1 at the correponding index
  if (l == -1)
	mutations = rand(N, var_count, 1) <= Pm;
  else
	mutations = rand(N, l, var_count) <= Pm;
  end

  %% Mutation
  result = mutation_fn(population, mutations, context);
end

function result = make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context)
  children = crossover_(mating_pool, crossover_fn, Pc, context);
  result = mutate_(children, l, mutation_fn, Pm, context);
end


function [all_h, plot_legend] = iterativeShowPaleto_(problem, variables, plot_optimal, algo_name, all_h, plot_legend)
  global UTILS;
  
  if (~exist('plot_optimal', 'var'))
      plot_optimal = false;
  end

  if (~exist('algo_name', 'var'))
      algo_name = 'Pareto Front';
  end
  
  if (~exist('all_h', 'var'))
      all_h = [];
  end
  
  if (~exist('plot_legend', 'var'))
      plot_legend = cell(0);
  end
  
  colors = ['r', 'b', 'g', 'k'];
  shape = ['*', 'd', '+', '.'];
  
  objective_vector = problem.objective_vector;

  [~, fn_count] = size(objective_vector); 

  can_plot_fn = (fn_count <= 3);
  
  if (~can_plot_fn)
	warning('showPaleto: variables: can not plot in more than 3D, sorry!');
  else
    hold on;
    
    if (plot_optimal && isempty(all_h))    
        pareto_front = UTILS.evalFnVector(objective_vector, problem.optimal_solutions(1000));
        
        if (fn_count == 3)
            plot_fn = @plot3;
        else
            plot_fn = @plot;
        end
        
        all_h(end+1) = plotN_(pareto_front, plot_fn, 'k-', 'LineWidth', 1);
        plot_legend{end+1} = 'Optimal pareto front';
    end
    
    values = UTILS.evalFnVector(objective_vector, variables);
       
    if (fn_count == 3)
        plot_fn = @plot3;
        zlabel("f3(X)");
    else
        plot_fn = @plot;
    end
  
    xlabel("f1(X)");
    ylabel("f2(X)");
    title("Objective's domain");
    
    i = length(all_h) + ~plot_optimal;
    style = sprintf('%c%c', colors(i), shape(i));
    
    all_h(end+1) = plotN_(values, plot_fn, style);
    
    plot_legend{end+1} = algo_name;
    
    legend(all_h, plot_legend);
    title(sprintf('Pareto front on %s', problem.name));
  end
end

function showPaleto_(problem, variables, plot_optimal, algo_name)
  global UTILS;
  figure(42)
  clf;
  hold on;
  
  if (~exist('plot_optimal', 'var'))
      plot_optimal = false;
  end

  if (~exist('algo_name', 'var'))
      algo_name = 'Pareto Front';
  end
  
  all_h = [];
  plot_legend = cell(0);
  
  objective_vector = problem.objective_vector;

  [~, fn_count] = size(objective_vector); 

  can_plot_fn = (fn_count <= 3);
  
  if (~can_plot_fn)
	warning('showPaleto: variables: can not plot in more than 3D, sorry!');
  else
    hold on;
    
    if (plot_optimal)    
        pareto_front = UTILS.evalFnVector(objective_vector, problem.optimal_solutions(1000));
        
        if (fn_count == 3)
            plot_fn = @plot3;
        else
            plot_fn = @plot;
        end
        
        all_h(end+1) = plotN_(pareto_front, plot_fn, 'k-', 'LineWidth', 1);
        plot_legend{end+1} = 'Optimal pareto front';
    end
    
    values = UTILS.evalFnVector(objective_vector, variables);
       
    if (fn_count == 3)
        plot_fn = @plot3;
        zlabel("f3(X)");
    else
        plot_fn = @plot;
    end
  
    xlabel("f1(X)");
    ylabel("f2(X)");
    title("Objective's domain");
    
    style = 'r*';
    all_h(end+1) = plotN_(values, plot_fn, style);
    
    plot_legend{end+1} = algo_name;
    
    legend(all_h, plot_legend);
    title(sprintf('Pareto front on %s', problem.name));
  end
end

function result = plotN_(val_array, plot_fn, varargin)
    BY_COLUMN = 2; 
    to_var_arg = num2cell(val_array', BY_COLUMN);
    
    result = plot_fn(to_var_arg{:}, varargin{:});
end

%% NOTE: This is needed to allow evaluation of objective_vector on
%% matrices (which is not handle with UTILS.evalFn (which was not
%% designed to handle it anyway).
%% Btw, this is an ugly hack!)
function result = evalFn_(objective_vector, varargin)
  [~, fn_count] = size(objective_vector);

  %% It works. I don't really know why, by it does!
  [N, M] = size(varargin{1});
  result = zeros(M, N, fn_count);

  for i = 1:fn_count
	f = objective_vector{i};

	result(:, :, i) = f(varargin{:});
  end
end

