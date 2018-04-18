function IBEA_ADAPTIVE
  global IBEA_ADAPTIVE;

  IBEA_ADAPTIVE.defaultConfig = @defaultConfig_;
  IBEA_ADAPTIVE.run = @run_;
  IBEA_ADAPTIVE.name = "IBEA_ADAPTIVE";
end

function [result, h] = run_(current_population, ga_context, config)
  global SELECTION;
  global GA;
  
  N = config.N;
  l = config.l;
  
  kappa = config.kappa;
  
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
    scale_objective_values = calculateScaleObjectiveValues(current_population, objective_vector, decode_fn);
    [fitness,maxAbsIndicatorValue] = calculateFitness(scale_objective_values, kappa);
    
    %step 3
    [sizePop,~] = size(current_population);    
    while(sizePop > N)
        
        posElementToRemove = findWorstIndex(fitness);
        [fitness, current_population, scale_objective_values] = removeAndRevaluate(fitness, current_population, scale_objective_values, posElementToRemove,maxAbsIndicatorValue, kappa);
        
        [sizePop,~] = size(current_population);
    end
    
    %step 4
    done = ((g == G_max)|| stop_criteria_fn(scale_objective_values, old_objective_values, maximizing));
    
    %step 5
    if (done)
        result = decode_fn(current_population);
        %save history
        h(g).population = result;
        h(g).objective_values = scale_objective_values;
    else
        selection = selection_fn(fitness);
        mating_pool = current_population(selection, :);
        
        new_population = GA.make_new_pop(mating_pool, l, crossover_fn, Pc, mutation_fn, Pm, context);
        new_population = GA.make_new_pop(current_population, l, crossover_fn, Pc, mutation_fn, Pm, context);
        current_population = [new_population ; current_population] ;% concatenate both populations
    
        old_objective_values = scale_objective_values;
        
        %save history
        h(g).population = current_population;
        h(g).objective_values = scale_objective_values;
        
        g = g + 1;
    end
  end
end


function scale_objective_values = calculateScaleObjectiveValues(population, fn_vector, decode_fn)
    global UTILS;
    
    real_values_pop = decode_fn(population);
    objective_values = UTILS.evalFnVector(fn_vector, real_values_pop);
    [N, fn_count] = size(objective_values);

    upperBound = zeros(fn_count);
    lowerBound = inf(fn_count);
    %find lower and upper bound
    for i = 1:fn_count
        for j = 1:N
            upperBound(i) = max(upperBound(i), objective_values(j,i) );
            lowerBound(i) = min(lowerBound(i), objective_values(j,i) );
        end
    end
    %create a scale objective value ( scale_objective_values(j,i) =< 1 for all i,j )
    scale_objective_values = zeros(N, fn_count);
    for i = 1:fn_count
        for j = 1:N
            scale_objective_values(j,i) = ( (objective_values(j,i) - lowerBound(i)) / ( upperBound(i) - lowerBound(i) ) );
        end
    end
end

function [fitness, maxAbsIndicatorValue] = calculateFitness(scale_objective_values, kappa)
    
    [N, fn_count] = size(scale_objective_values);
    fitcomp = zeros(N,N);
    maxAbsIndicatorValue = -inf;
    
    for i = 1:N
        for j = 1:N
            fitcomp(i,j) = additiveIndicator(scale_objective_values(i,:), scale_objective_values(j,:) , fn_count);  
            if (abs(fitcomp(i,j)) > maxAbsIndicatorValue)
                maxAbsIndicatorValue = abs( fitcomp(i,j) );
            end
        end
    end
    
    for i = 1:N
        sum = 0;
        for j = 1:N
            if ( i ~= j)
                sum = sum +  exp((-fitcomp(i,j) / maxAbsIndicatorValue ) / kappa);
            end
        end
        fitness(i) = sum;
    end
end

%
function [fitness, current_population, scale_objective_values] = removeAndRevaluate(fitness, current_population, scale_objective_values, posElementToRemove,maxAbsIndicatorValue, kappa)
    [~, fn_count] = size(scale_objective_values);
    [sizePop,~] = size(current_population);

    for i = 1:sizePop
        if (i ~= posElementToRemove)
            epsilon = additiveIndicator(scale_objective_values(i,:), scale_objective_values(posElementToRemove,:) , fn_count);     
            fitness(i) = fitness(i) - exp((-epsilon/maxAbsIndicatorValue ) / kappa);
        end
    end
    
    current_population(posElementToRemove,:) =[]; %remove this element
    fitness(posElementToRemove) =[]; %remove this element
    scale_objective_values(posElementToRemove, :) =[]; %remove this element
end

function epsilon = additiveIndicator(a, b, fn_count)
    epsilon = -inf;
    for i = 1:fn_count
        epsilon = max ( a(i)-b(i), epsilon );
    end
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