function StopCriteria
			%STOPCRITERIA All stop criteria functions.
			%
			% time
			% threshold(RELATION, T), T as limit (depends on RELATION)
			% variance(V), V as lower limit
			% minMaxRatio(R)
			% meanChangeRate(CR), CR as lower limit
			%
			% See also STOPCRITERIA>THRESHOLD, STOPCRITERIA>VARIANCE, 
			% STOPCRITERIA>MINMAXRATIO, STOPCRITERIA>MEANCHANGERATE
  global STOP_CRITERIA;
  
  %% Time
  STOP_CRITERIA.time = @(f, old_f) 0;
  
  STOP_CRITERIA.threshold = @threshold;
  STOP_CRITERIA.variance = @variance;
  STOP_CRITERIA.minMaxRatio = @minMaxRatio;

  STOP_CRITERIA.meanChangeRate = @meanChangeRate;
end

function h = threshold(relation, t)
			%THRESHOLD Stop at least one individual as a fitness such as RELATION(FITNESS, T) == 1.
  
  h = @(f, old_f) threshold_(relation, t, f);
end

function result = threshold_(relation, t, fitness)
  result = ~isempty(find(relation(fitness, t), 1));
end

function h = variance(v)
  %VARIANCE Stop when var(fitness) <= V..
  
  h = @(f, old_f) variance_(v, f);
end

function result = variance_(v, fitness)
  result = (var(fitness) <= v);
end

function h = minMaxRatio(r)
 %MINMAXRATIO Stop when max(fitness) / min(fitness) is equal to ratio.
  
  if (r <= 0)
	error('minMaxRatio: R must be > 0');
  end
  
  h = @(f, old_f) minMaxRatio_(r, f);
end

function result = minMaxRatio_(ratio, fitness)
  max_f = max(fitness);
  min_f = min(fitness);

  result = (abs((max_f / min_f) - ratio) <= eps);
end

function h = meanChangeRate(cr)
	%MEANCHANGERATE Stop when the difference between mean(fitness) and
	% mean(old_fitness) is <= CR * 100.
  
  h = @(f, old_f) meanChangeRate_(cr, f, old_f);
end

function result = meanChangeRate_(change_rate, fitness, old_fitness)
  %% First iteration, no old fitness
  if (length(old_fitness) == 0)
	result = 0;
	return;
  end
  
  mean_f = mean(fitness);
  mean_old_f = mean(old_fitness);

  rate = abs((mean_f - mean_old_f) / mean_old_f);
  result = (rate <= change_rate);
end
