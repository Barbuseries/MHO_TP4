function StopCriteria
			%STOPCRITERIA All stop criteria functions.
			%
			% time
			% threshold(T), T as limit (upper limit if minimizing, upper limit otherwhise)
			% variance(V), V as lower limit
			% minMaxRatio(R)
			% meanChangeRate(CR), CR as lower limit
			%
			% See also STOPCRITERIA>THRESHOLD, STOPCRITERIA>VARIANCE, 
			% STOPCRITERIA>MINMAXRATIO, STOPCRITERIA>MEANCHANGERATE
  global STOP_CRITERIA;
  
  %% Time
  STOP_CRITERIA.time = @(varargin) 0;
  
  STOP_CRITERIA.threshold = @threshold;
  STOP_CRITERIA.variance = @variance;
  STOP_CRITERIA.minMaxRatio = @minMaxRatio;

  STOP_CRITERIA.meanChangeRate = @meanChangeRate;
end

function h = threshold(t, count)
%THRESHOLD Stop at least COUNT individuals have all objective values <= T
% if minimizing, >= T otherwhise.
%  
% T can either be a scalar (same T for all objective values) or a
% vector.
%
% COUNT defaults to 1.

  if (~exist('count', 'var'))
    count = 1;
  elseif (count <= 0)
    error('threshold: COUNT must be greater than 0');
  end
  
  h = @(ov, ~, maximizing) threshold_(t, count, ov, maximizing);
end

function result = threshold_(t, count, objective_values, maximizing)
  [N, fn_count] = size(objective_values);
  t_count = length(t);
  
  if (t_count == 1)
	t = t * ones(1, fn_count);
  else
	if (t_count ~= fn_count)
	  error('threshold: T must either be a scalar or a vector the same length as the OBJECTIVE_VECTOR');
	end
  end
  
  if (count > N)
      error('threshold: COUNT must be less than or equal to N');
  end
  
  if (maximizing)
	relation = @ge;
  else
	relation = @le;
  end

  BY_ROW = 2;
  is_valid = sum(relation(objective_values, t), BY_ROW) == fn_count;
  result = sum(is_valid) >= count;
end

function h = variance(v)
  %VARIANCE Stop when var(objective_values) <= V..
  
  h = @(ov, varargin) variance_(v, ov);
end

%% TODO: Same as t, allow scalar or vector.
function result = variance_(v, objective_values)
  [~, fn_count] = size(objective_values);
  result = sum(var(objective_values) <= v) == fn_count;
end

function h = minMaxRatio(r)
%MINMAXRATIO Stop when max(objective_values) / min(objective_values)
% is equal to ratio.
  
  if (r <= 0)
	error('minMaxRatio: R must be > 0');
  end
  
  h = @(ov, varargin) minMaxRatio_(r, ov);
end

%% TODO: Same as t, allow scalar or vector.
function result = minMaxRatio_(ratio, objective_values)
  [~, fn_count] = size(objective_values);
  
  max_ov = max(objective_values);
  min_ov = min(objective_values);

  result = sum(abs((max_ov ./ min_ov) - ratio) <= eps) == fn_count;
end

%% TODO: Same as t, allow scalar or vector.
function h = meanChangeRate(cr)
%MEANCHANGERATE Stop when the difference between mean(objective_values) and
% mean(old_objective_values) is <= CR * 100.
  
  h = @(ov, old_ov, ~) meanChangeRate_(cr, ov, old_ov);
end

function result = meanChangeRate_(change_rate, objective_values, old_objective_values)
  [~, fn_count] = size(old_objective_values);

  %% First iteration, no old objective_values
  if (fn_count == 0)
	result = 0;
	return;
  end
  
  mean_ov = mean(objective_values);
  mean_old_ov = mean(old_objective_values);

  rate = abs((mean_ov - mean_old_ov) ./ mean_old_ov);
  result = sum(rate <= change_rate) == fn_count;
end
