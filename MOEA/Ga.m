function Ga
  global GA;
  GA.showPaleto = @showPaleto_;
  GA.defaultConfig = @defaultConfig;
  GA.make_new_pop = @make_new_pop;
  
  GA.create = @create_ga_;
end

function result = optimize(ga, maximizing, objective_vector, constraints, config)
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
  
  tic;
  
  population = initialGeneration_(N, constraints, l);
  result = ga.run(population, ga_context, config);
  
  toc;
end

function result = maximize(ga, objective_vector, constraints, config)
  result = optimize(ga, 1, objective_vector, constraints, config);
end

function result = minimize(ga, objective_vector, constraints, config)
  result = optimize(ga, 0, objective_vector, constraints, config);
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

function result = create_ga_(ga)
  result = ga;
  
  result.optimize = @(varargin) optimize(ga, varargin{:});
  result.maximize = @(varargin) maximize(ga, varargin{:});
  result.minimize = @(varargin) minimize(ga, varargin{:});
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

function result = make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context)
  [N, var_count] = size(mating_pool);
  
  children = crossover_(mating_pool, crossover_fn, Pc, context);

  %% Every allele that needs to mutate is 1 at the correponding index
  if (l == -1)
	mutations = rand(N, var_count, 1) <= Pm;
  else
	mutations = rand(N, l, var_count) <= Pm;
  end

  %% Mutation
  result = mutation_fn(children, mutations, context);
end


%% TODO: Test 3D function plot
function showPaleto_(problem, variables, plot_objective_domain)
  global UTILS;
  figure(5)
  clf;
  
  if (~exist('plot_objective_domain', 'var'))
      plot_objective_domain = false;
  end

  constraints = problem.constraints;
  objective_vector = problem.objective_vector;

  [~, var_count] = size(variables);
  [~, fn_count] = size(objective_vector);

  can_plot_var = (var_count <= 3);
  can_plot_fn = (fn_count <= 3);
  plot_count = can_plot_var + can_plot_fn;
  
  if (var_count == 1)
      m = 500;
  elseif (var_count == 2)
      m = 50;
  else
      m = 20;
  end
  
  %% FIXME: This does not work if some constraints are different from the others
  if (var_count == 1)
    dx = meshgrid(UTILS.linspacea(constraints, m));  
    dy = 0 * dx;
    dz = 0 * dy;
  elseif (var_count == 2)
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

    plotted_values = variables;
    plot_axis = reshape(constraints', 1, []);
    
    %% Set plot's axis limits to constraints
    if (var_count == 1)    
        %% Set a fictional y to 0, for a prettier plot
        plotted_values(:, end+1) = 0;
        plot_axis = [plot_axis, -0.1, 0.1];
    end
    
    axis(plot_axis);
    
    if (var_count == 3)
        plot_fn = @plot3;
        zlabel("z");
    else    
        plot_fn = @plot;
    end
    
    xlabel("x");
    ylabel("y");
    title("Variables' space");
    plotN_(plotted_values, plot_fn, 'r+');
  end

  if (~can_plot_fn)
	warning('showPaleto: variables: can not plot in 4D, sorry!');
  else
	subplot(1, plot_count, plot_count);
    hold on;

    if (plot_objective_domain)
        if (var_count == 1)
            fz = evalFn_(objective_vector, dx);
        elseif (var_count == 2)
            fz = evalFn_(objective_vector, dx, dy);
        else
            fz = evalFn_(objective_vector, dx, dy, dz);
        end

        %% TODO: Test domain plot in 3D (no 3D function yet)
        if (fn_count == 3)
            mesh(fz(:, :, 1), fz(:, :, 2), fz(:, :, 3));
        else
            x = fz(:, :, 1);
            y = fz(:, :, 2);
            z = zeros(size(x));
            C = gradient(x .* y);

            mesh(x, y, z, C);
        end
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
    h = plotN_(values, plot_fn, 'r+');
    
    legend(h, "Pareto front");
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

