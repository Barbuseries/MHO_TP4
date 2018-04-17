function IBEA
  global IBEA;

  IBEA.defaultConfig = @defaultConfig_;
  IBEA.run = @run_;
  IBEA.name = "IBEA";
end

function [result, h] = run_(current_population, ga_context, config)
  global SELECTION;
  global GA;
  global UTILS;
  
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
  g = 1;

  h_template = struct('population', [], 'objective_values', []);
  h = repmat(h_template, G_max, 1);
  
  old_objective_values = [];
  archive = [];
  done = false;
  
  kappa = 1; %(fitness scaling factor)
  
  while (~done)
	if (l == -1)
	  context.iteration = g;
    end
    
    %step 2
    [fitness, real_values_pop, objective_values, maxAbsIndicatorValue] = evalFitness_(current_population, objective_vector, decode_fn, kappa);
    
    %step 3
    [sizePop,~] = size(current_population);    
    
    while(sizePop > N)
        posElementToRemove = findWorstIndex(fitness);
        [fitness, current_population, objective_values] = removeAndRevaluate(fitness, current_population, objective_values, posElementToRemove, maxAbsIndicatorValue, kappa);
        
        [sizePop,~] = size(current_population);
    end
    
    %step 4
    done = ((g == G_max)|| stop_criteria_fn(objective_values, old_objective_values, maximizing));
    
    %step 5
    if (done)
        result = real_values_pop;
        %save history
        h(g).population = result;
        h(g).objective_values = objective_values;
    else
        selection = selection_fn(fitness);
        mating_pool = current_population(selection, :);
        
        new_population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);
        current_population = [new_population ; current_population] ;% concatenate both populations
    
        old_objective_values = objective_values;
        
        %save history
        h(g).population = current_population;
        h(g).objective_values = objective_values;
        
        g = g + 1;

    end
  end
end

function [fitness, real_values_pop, objective_values, maxAbsIndicatorValue] = evalFitness_(population, fn_vector, decode_fn, kappa)
    global UTILS;
    
    real_values_pop = decode_fn(population);
    objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);
    [N, fn_count] = size(objective_values);
    
    %find maxAbsIndicatorValue
    maxAbsIndicatorValue = -inf;
    fitcomp = zeros(N,N);
    for i = 1:N
        for j = 1:N
            fitcomp(i,j) = additiveIndicator(objective_values(i,:), objective_values(j,:) , fn_count);
            if ( abs(fitcomp(i,j)) > maxAbsIndicatorValue)
                maxAbsIndicatorValue = abs(fitcomp(i,j));
            end
            
        end
    end
    %calculate fitness
    fitness = zeros(1,N);
    for i = 1:N
        sum = 0;
        for j = 1:N
            if (i ~= j)
                sum = sum + exp((-fitcomp(i,j)/maxAbsIndicatorValue ) / kappa);
            end
        end
        fitness(i) = sum;
    end
end

function epsilon = additiveIndicator(a, b, fn_count)
    epsilon = -inf;
    for i = 1:fn_count
        epsilon = max ( a(i)-b(i), epsilon );
    end
end

function [fitness, current_population, objective_values] = removeAndRevaluate(fitness, current_population, objective_values, posElementToRemove, maxAbsIndicatorValue, kappa)
    [~, fn_count] = size(objective_values);
    [sizePop,~] = size(current_population);

    for i = 1:sizePop
        if (i ~= posElementToRemove)
            epsilon = additiveIndicator(objective_values(i,:), objective_values(posElementToRemove,:) , fn_count);     
            fitness(i) = fitness(i) - exp((-epsilon/ maxAbsIndicatorValue ) / kappa);
        end
    end
    current_population(posElementToRemove,:) =[]; %remove this element
    fitness(posElementToRemove) =[]; %remove this element
    objective_values(posElementToRemove, :) =[]; %remove this element
end

function result = findWorstIndex(fitness)
    %create ranking
    [~,~,rank] = unique(fitness);
    rank = rank';
    %select the ind with the smalless fitness
    result = find( rank == 1, 1 );
end

function result = defaultConfig_
  global GA;
  result = GA.defaultConfig();
  
end


