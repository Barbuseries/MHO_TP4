function Ibea
  global IBEA;

  IBEA.defaultConfig = @defaultConfig_;
  IBEA.run = @run_;
  IBEA.name = "IBEA";
end

function [result, h] = run_(current_population, ga_context, config)
  global SELECTION;
  global GA;
  
  N = config.N;
  l = config.l;
  
  kappa = config.kappa;
  adaptative = config.adaptative;
  
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
  g = 1;

  h_template = struct('population', [], 'objective_values', []);
  h = repmat(h_template, G_max, 1);
  
  old_objective_values = [];
  done = false;
  
  while (~done)
	if (l == -1)
	  context.iteration = g;
    end

    %step 2
    [fitness, ind_exp, real_values_pop, objective_values] = evalFitness_(current_population, objective_vector, decode_fn, kappa, adaptative);
    
    %step 3
    if (g > 1)
	  kept_indices = 1:(2*N);
    
	  for i = 1:N
		[~, posElementToRemove] = min(fitness(kept_indices));
        j = kept_indices(posElementToRemove);
        
        kept_indices(posElementToRemove) = [];
        
		fitness = removeAndRevaluateIterative(fitness, posElementToRemove, ind_exp, kept_indices);
      end

	  fitness = fitness(kept_indices);
      current_population = current_population(kept_indices, :);
      objective_values = objective_values(kept_indices, :);
    end
    
    %step 4
    done = ((g == G_max)|| stop_criteria_fn(objective_values, old_objective_values, maximizing));

	non_dominated = get_non_dominated_(objective_values, maximizing);
	
	%%save history
    h(g).population = current_population(non_dominated, :);
    h(g).objective_values = objective_values(non_dominated, :);
    
    %step 5
    if (done)
        %% Even though those elements were removed from the
        %% population, we never removed them from this...
        real_values_pop = real_values_pop(kept_indices, :);
        
        result = real_values_pop(non_dominated, :);
        associated_values = objective_values(non_dominated, :);
        
        %save history
        h(g).population = result;
        h(g).objective_values = associated_values;
    else
        %% NOTE: It _seems_ fitness is better when >.
        selection = selection_fn(fitness);
        mating_pool = current_population(selection, :);
        
        new_population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);
        current_population = [new_population ; current_population] ;% concatenate both populations
    
        old_objective_values = objective_values;
        
        g = g + 1;
    end
  end
end

function [fitness, ind_exp, real_values_pop, objective_values] = evalFitness_(population, fn_vector, decode_fn, kappa, adaptative)
    global UTILS;

	BY_COLUMN = 1;
    BY_ROW = 2;
    
    real_values_pop = decode_fn(population);
    objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);
    [N, ~] = size(objective_values);

	%% Normalize objective values
	if (adaptative)
	  min_bounds = min(objective_values, [], BY_COLUMN);
	  max_bounds = max(objective_values, [], BY_COLUMN);

	  values = (objective_values - min_bounds) ./ (max_bounds - min_bounds);
	else
	  values = objective_values;
	end
	
    %% For each individual i, compute epsilon(j, i), which is the
    %% minimum amount by which a plane, parallel to one of the axes
    %% and crossing j, needs to move to exceed or equal i. (epsilon is
    %% the maximum amount between i and j in any coordinate).
    epsilon = zeros(N, N);
    for i = 1:N
      epsilon(i, :) = additiveIndicator(values, values(i,:));
    end

	%% Normalize epsilon
	if (adaptative)
	  epsilon = epsilon / max(abs(epsilon(:)));
	end
	
    %calculate fitness
    ind_exp = exp(-epsilon / kappa);
	
	%% +1, because we always compute epsilon(i, i), which gives us
	%% -e^(0) (-1) in the end. This is 100% not useful (because all
	%% fitnesses are offset by it anyway...).
    fitness = (-sum(ind_exp, BY_ROW) + 1)';
end

function epsilon = additiveIndicator(a, b)
  %% @perf Btw, this is a one-liner, and it is the bottleneck. Houray!
  
    BY_ROW = 2;
    epsilon = max(a - b, [], BY_ROW)';
end

function [fitness] = removeAndRevaluateIterative(fitness, posElementToRemove, ind_exp, kept_indices)
    BY_ROW =  2;
        
    %% For each removed element, get back the associated indicator
    %% value for each kept individual.
    %% Add this value to the individual's fitness.
    %% (NOTE: We already store the exp^(-I(x, i) / kappa) row-wise,
    %% meaning that we have a matrix where every row is exp^(-I(x, i)
    %% / kappa).
    %% This is tricky: we need exp^(-I(x*, x) / kappa), where x* is
    %% the removed individual, so we first get the rows we care about
    %% (individuals that are not removed, all 'x's), then the
    %% columns. As each row corresponds to a given individual, we just
    %% sum the values to get the it's actual delta fitness.)
	
    relevent_indicators = ind_exp(kept_indices, posElementToRemove);
    fitness(kept_indices) = fitness(kept_indices) + sum(relevent_indicators, BY_ROW)';
end

function non_dominated = get_non_dominated_(objective_values, maximizing)
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
	vals = objective_values(i, :);
	dominators = (sum(relation(objective_values, vals), BY_ROW) == fn_count) & (sum(objective_values == vals, BY_ROW) < fn_count);
    
	result(i) = sum(dominators) ~= 0;
  end
end


function result = defaultConfig_
  global GA;
  result = GA.defaultConfig();

  result.kappa = 1; %% fitness scaling factor
  result.adaptative = false;
end
