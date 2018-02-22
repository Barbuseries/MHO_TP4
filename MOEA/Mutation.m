function Mutation
	   %MUTATION All mutation functions.
	   %
	   % Binary mutations
	   %  bitFlip
	   % 
	   % Arithmetic mutations
	   %  uniform
	   %  boundary
	   %  normal(SIGMA), length(SIGMA) is either 1 (same SIGMA for all
	   %    variables) or VAR_COUNT
	   %  polynomial(N), N >= 0
	   %  nonUniform(B)
	   %
	   % See also MUTATION>BITFLIP, MUTATION>UNIFORM,
	   % MUTATION>BOUNDARY, MUTATION>NORMAL, MUTATION>POLYNOMIAL,
	   % MUTATION>NONUNIFORM
  global MUTATION;

  %% Binary
  MUTATION.bitFlip = @bitFlip;

  %% Arithmetic
  MUTATION.uniform = @uniform;
  MUTATION.boundary = @boundary;
  MUTATION.normal = @normal;
  MUTATION.polynomial = @polynomial;
  MUTATION.nonUniform = @nonUniform;
end

%% Binary
function result = bitFlip(children, mutations, ~)
	 %BITFLIP For each element I in mutations where mutations(I) == 1,
	 % invert the corresponding bit in CHILDREN.
  
  global UTILS;
  
  dim = size(children);
  var_count = dim(2);

  %% 'mutations' contains 1s and 0s, that we use to represent the
  %% mutation (1s flip the bit when xor-ing, 0s do not).
  %% 'reshape' is used to handle different masks for each variable.
  %% (We have a matrix of Nxlxvar_count 1s and 0s, that we them
  %% convert to Nxvar_count integers)
  mask = reshape(UTILS.arrayToDec(mutations), [], var_count);

  result = bitxor(children, mask); %% Do a flip!
end

%% Arithmetic
function result = uniform(children, mutations, context)
	 %UNIFORM For each element I in mutations where mutations(I) == 1,
	 % assign a random value inside CONTEXT.CONSTRAINTS to the
	 % corresponding variable in children.
  
  global UTILS;
    
  dim = size(children);

  N = dim(1);
  
  %% Children which do not mutate (mutations == 0) keep their values.
  %% Otherwhise, they get a random value inside the constraints.
  result = children .* (mutations == 0) + mutations .* UTILS.randomIn(context.constraints, N);
end

function result = boundary(children, mutations, context)
	%BOUNDARY For each element I in mutations where mutations(I) == 1,
	% draw a random number U in ]0, 1[. If assign U > 0.5, set the
	% corresponding variable in children to the upper bound in
	% CONTEXT.CONSTRAINTS, otherwhise, set it to the lower bound.
  
  constraints = context.constraints;
  
  dim = size(children);

  N = dim(1);
  var_count = dim(2);
  
  %% If we do not mutate, mutation is 0.
  %% Otherwhise, it is set to a random number in ]0, 1[.
  mutations = mutations .* rand(N, var_count, 1); 
  
  keep = children .* (mutations == 0); %% Keep value (or zero if not related)
  clamp_up = constraints(:, 2)' .* (mutations > 0.5); %% Set to upper bound (or zero if not related)
  clamp_down = constraints(:, 1)' .* ((mutations > 0) & (mutations <= 0.5)); %% Set to lower bound (or zero if not related)
  
  result = keep + clamp_up + clamp_down;
end

function h = normal(sigma)
  %NORMAL Return a function that produces
  % NORMAL_(SIGMA, CHILDREN, MUTATIONS, CONTEXT), when given
  % CHILDREN, MUTATIONS and CONTEXT.
  %   H = NORMAL(SIGMA)
  %
  % See also MUTATION>NORMAL_.
  h = @(c, m, cx) normal_(sigma, c, m, cx);
end

function result = normal_(sigma, children, mutations, context)
  %NORMAL_ Normal mutation. SIGMA can either be the same for all
  %variables (in tha case, a scalar is enough), or different for each
  %(an array with as many elements as there are variables is needed)
  %
  % See also MUTATION>NORMAL, NORMAN>NORMAL.
  constraints = context.constraints;
  clamp_fn = context.clamp_fn;
  [~, var_count] = size(children);

  sigma_len = length(sigma);

  if (sigma_len == 1)
	%% Use the same sigma for all variables
	sigma = sigma * ones(1, var_count);
  elseif (sigma_len ~= var_count)
	error('normal: length(SIGMA) must either be 1 or equal to VAR_COUNT');
  end

  %% Sigma is set to 0 if there is no mutation
  sigma = sigma .* (mutations == 1);

  non_zero = (sigma ~= 0);
  sigma_non_zero = sigma(non_zero);
  sigma(non_zero) = sigma_non_zero .* normrnd(0, 1, size(sigma_non_zero));
  
  lowest = constraints(:, 1)';
  biggest = constraints(:, 2)';
  
  result = children + sigma;
  result = clamp_fn(result, lowest, biggest);
end

function h = polynomial(n)
			%POLYNOMIAL Return a function that produces
			% POLYNOMIAL_(N, CHILDREN, MUTATIONS, CONTEXT), when given
			% CHILDREN, MUTATIONS and CONTEXT.
			%   H = POLYNOMIAL(N)
			%
			% See also MUTATION>POLYNOMIAL_.

  if (n < 0)
	error('polynomial: N must be >= 0');
  end
  
  h = @(c, m, cx) polynomial_(n, c, m, cx);
end


function result = polynomial_(n, children, mutations, context)
							%POLYNOMIAL_ Polynomials and stuff. It works.
							%
							% See also MUTATION>POLYNOMIAL.
  
  constraints = context.constraints;
  clamp_fn = context.clamp_fn;
  
  dim = size(children);
  
  N = dim(1);
  var_count = dim(2);
  
  delta_max = (constraints(:, 2) - constraints(:, 1))';

  non_zero = (mutations == 1);
  u = rand(N, var_count);

  u_below = (u < 0.5);
  u_above = ~u_below;

  %% Do the computation on u as a whole, but apply the condition after
  %% (by multiplying by either 1 or 0 if the condition was true or
  %% not).
  %% This is _way_ faster than a for-loop.
  inv = 1 / (n + 1);
  xi = ((2 * u).^inv  - 1) .* u_below + (1 - (2 * (1 - u)).^inv) .* u_above;

  lowest = constraints(:, 1)';
  biggest = constraints(:, 2)';
  
  result = children + delta_max .* xi .* non_zero;
  result = clamp_fn(result, lowest, biggest);
end

function h = nonUniform(b)
			%NONUNIFORM Return a function that produces
			% NONUNIFORM_(N, CHILDREN, MUTATIONS, CONTEXT), when given
			% CHILDREN, MUTATIONS and CONTEXT.
			%   H = NONUNIFORM(N)
			%
			% See also MUTATION>NONUNIFORM_.
  
  h = @(c, m, cx) nonUniform_(b, c, m, cx);
end

function result = nonUniform_(b, children, mutations, context)
								%NONUNIFORM_ It is non-uniform alright.
								%
								% See also MUTATION>NONUNIFORM.
  
  constraints = context.constraints;
  clamp_fn = context.clamp_fn;
  
  g = context.iteration;
  G_max = context.G_max;

  dim = size(children);
  N = dim(1);
  var_count = dim(2);

  lowest = constraints(:, 1)';
  biggest = constraints(:, 2)';

  non_zero = (mutations == 1);
  u = rand(N, var_count);

  u_below = (u < 0.5);
  u_above = ~u_below;

  %% Do the computation on u as a whole, but apply the condition after
  %% (by multiplying by either 1 or 0 if the condition was true or
  %% not).
  %% This is _way_ faster than a for-loop.
  inv = ((1 - g) / G_max)^b;
  delta_g = (1 - u.^inv);
  
  delta_above = (biggest - children) .* delta_g .* u_above;
  delta_below = (children - lowest) .* delta_g .* u_below;

  result = children + (delta_above - delta_below) .* non_zero;
  result = clamp_fn(result, lowest, biggest);
end
