function Benchmark
  global BENCHMARK;
  
  BENCHMARK.run = @run_;
  BENCHMARK.plot = @plot_;
end

function result = run_(all_problems, all_algos, all_configs, count)
  global GA;
  global UTILS;
  
  problem_count = length(all_problems);
  algo_count = length(all_algos);

  result_template = struct('metrics', ...
						   struct('convergence', [], ...
                                  'diversity', []));
  
  result = repmat(result_template, algo_count, problem_count, count);
  total_iteration_count = problem_count * count * algo_count;
  iteration = 1;
  start = datetime('now');

  for k = 1:problem_count
	problem = all_problems{k};
    p_dummy = problem([]);
	eval_problem = @(x) UTILS.evalFnVector(p_dummy.objective_vector, x);
	pareto_front = eval_problem(p_dummy.optimal_solutions(500));
	
	for i = 1:count
      initial_population = [];
      all_Ns = cellfun(@(c) c.N, all_configs(:, k));
      all_ls = cellfun(@(c) c.l, all_configs(:, k));
      
      all_equal = @(a) sum(a) == (a(1) * length(a));
      if (all_equal(all_Ns) && all_equal(all_ls))
          N = all_Ns(1); l = all_ls(1);
          initial_population = GA.initialPopulation(N, p_dummy.constraints, l);
      end          
      
	  for j = 1:algo_count
        ratio = (iteration / total_iteration_count);
        estimated_time = NaN;
        delta_time = time(between(start, datetime('now')));
        if (iteration > 0)
            estimated_time = delta_time * ((1 - ratio) / ratio);
        end
        fprintf(1, '\t\tIteration: %d / %d (%.02f%%) - Estimated time: %ss\n', ...
                    iteration, total_iteration_count, 100 * ratio, estimated_time);
          
		algo = GA.create(all_algos(j));
		
		p = problem(algo);
		
		config = algo.defaultConfig();
		config = mergeStruct(config, all_configs{j, k});
		config.population = initial_population;

		[~, h] = p.optimize(config);

		result(j, k, i).metrics = computeMetrics_({h.objective_values}, pareto_front);
		
        iteration = iteration + 1;
	  end
	end
  end
end

function result = computeMetrics_(h, pareto_front)
  BY_ROW = 2;
  BY_COLUMN = 1;
  
  N = length(h);

  result = struct('convergence', zeros(1, N), 'diversity', zeros(1, N));
  [~, sorted_pareto_indices] = sort(pareto_front(:, 1));
  sorted_pareto_front = pareto_front(sorted_pareto_indices, :);
  
  for i = 1:N
	obtained_pareto_front = h{i};
	
  	if (isempty(obtained_pareto_front))
      result.convergence(i) = NaN;
    else
      [M, ~] = size(obtained_pareto_front);
      
      %% Convergence metric
      all_min_distances_sq = zeros(M, 1);
      for n = 1:M
		all_min_distances_sq(n) = min(sum((obtained_pareto_front(n, :) - pareto_front).^2, BY_ROW));
      end

      result.convergence(i) = mean(sqrt(all_min_distances_sq));
      
      %% Diversity metric
      if (M < 2)
        result.diversity(i) = NaN;
      else
        distance_fn = @(a, b) sqrt(sum(((b - a)).^2, BY_ROW));
          
        [~, sorted_indices] = sort(obtained_pareto_front(:, 1));
        sorted_values = obtained_pareto_front(sorted_indices, :);
        consecutive_distances = distance_fn(sorted_values(1:end-1, :), sorted_values(2:end, :));
        mean_consecutive_distance = mean(consecutive_distances);
      
        %% FIXME: This should values on a parallel curve instead of the pareto front.
        %% But I do not know how to compute that =(.
        d_f = distance_fn(sorted_values(1, :), sorted_pareto_front(1, :));
        d_l = distance_fn(sorted_values(end, :), sorted_pareto_front(end, :));
        
        upper = (d_f + d_l + sum(abs(consecutive_distances - mean_consecutive_distance)));
        lower = (d_f + d_l + ((M - 1) * mean_consecutive_distance));
      
        result.diversity(i) =  upper / lower;
      end
	end
  end
end

function plot_(all_algos, all_problems, data, compare_exact) 
  if (~exist('compare_exact', 'var'))
    compare_exact = false;
  end
  
  [algo_count, problem_count, run_count] = size(data);

  colors = ['r', 'b', 'g', 'k'];
  shape = ['*', 'd', '+', '.'];
  
  fields = fieldnames(data(1).metrics);
  field_count = length(fields);
	
  for j = 1:problem_count
	figure(j + 30);
	clf;
	hold on;
	
    p = all_problems{j}([]);

	for only_last_gen = 0:1
	  for k = 1:field_count
		component = fields{k};
		subplot(field_count, 2, field_count * only_last_gen + k);
		hold on;
		
		all_h = [];
		
		if (only_last_gen)
          a = [data(:, j, :).metrics];
          formated_data = zeros(run_count, algo_count);
          temp_data = reshape(a, algo_count, []);
          
          %% NOTE: I could probably one-line or two-line this, but it is starting to make my head hurt, so not today.
          for i = 1:algo_count
			temp_data2 = [temp_data(i, :).(component)];
			N = numel(temp_data2) / run_count;
			formated_data(:, i) = log10(temp_data2(N:N:end));
          end
%formated_data = (reshape([a.(component)], [], run_count, algo_count));
%formated_data = permute(formated_data(end, :, :), [3, 2, 1])';
%formated_data(isinf(formated_data)) = 0;
%formated_data(:, 1:end) = formated_data(:, end:-1:1);
          boxplot(formated_data, 'colors', colors, ...,
                  'Labels', cellstr([all_algos.name]), ...
                  'Notch', 'on');
          
		else
          for i = 1:algo_count
			a = [data(i, j, :).metrics];
			formated_data = log10(reshape([a.(component)], [], run_count)'); %% Reshape to have corresponding iterations on the same column.
							 %formated_data(isinf(formated_data)) = 0;

			[~, N] = size(formated_data);

			style = sprintf('%c%c%s', colors(i), shape(i), '-');

			if (compare_exact)
              boxplot(formated_data, 'colors', colors(i), 'Symbol', shape(i), ...
                      'PlotStyle', 'Compact', 'BoxStyle', 'outline');

              labels = 1:20:N;
            
              set(gca, 'XTick', labels);
              set(gca, 'XTickLabel', labels);
			else
              h = plot(1:N, mean(formated_data, 'omitnan'), style, 'MarkerSize', 4);
			end

			xlabel('Iteration');
			ylabel('Metric (Log10)');

			if (~compare_exact)
              all_h(end+1) = h;
			end
          end

          if (~compare_exact)
			legend(all_h, [all_algos.name]);
          end
		end

		if (only_last_gen)
		  title_format = 'Comparison of the %s metric on %s (result)';
		else
		  title_format = 'Comparison of the %s metric on %s';
		end

		title(sprintf(title_format, component, p.name));
      end
	end
  end
end
