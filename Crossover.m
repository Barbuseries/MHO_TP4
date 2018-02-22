function Crossover
	   %CROSSOVER All crossover functions.
	   %
	   % Binary crossovers
	   %  singlePoint
	   %  multiPoint(N), N in [1, L - 1]
	   %  uniform(P, T), P in [0, 1];  or P and T as function handles
	   %  oneBitAdaptation(F0, F1), F0 and F1 as crossover functions
	   % 
	   % Arithmetic crossovers
	   %  wholeArithmetic
	   %  localArithmetic
	   %  blend(ALPHA = 0.5)
	   %  simulatedBinary(N), N >= 0
	   %
	   % See also CROSSOVER>SINGLEPOINT, CROSSOVER>MULTIPOINT,
	   % CROSSOVER>ONEBITADAPTATION, CROSSOVER>UNIFORM,
	   % CROSSOVER>WHOLEARITHMETIC, CROSSOVER>LOCALARITHMETIC,
	   % CROSSOVER>BLEND, CROSSOVER>SIMULATEDBINARY
  
  global CROSSOVER;

  %% Binary
  CROSSOVER.singlePoint = @singlePoint;
  CROSSOVER.multiPoint = @multiPoint;
  CROSSOVER.uniform = @uniform;
  CROSSOVER.oneBitAdaptation = @oneBitAdaptation;


  %% Arithmetic
  CROSSOVER.wholeArithmetic = @wholeArithmetic;
  CROSSOVER.localArithmetic = @localArithmetic;
  CROSSOVER.blend = @blend;
  CROSSOVER.simulatedBinary = @simulatedBinary;
end

%% Binary crossovers
function result = singlePoint(a, b, l)
	%SINGLEPOINT Choose a random point in [1, L - 1] and split both A
	% and B into two parts, then swap and merge them to create two
	% children.
	% A different point is choosen for each variable in A and B.
	%   CHILDREN = SINGLEPOINT(A, B, L) Create two children from A and
	%   B, with a chromosome length of L.
	%
	% See also CROSSOVER>MULTIPOINT, CROSSOVER>MULTIPOINT_.
  
  dim = size(a);
  N = dim(1);
  var_count = dim(2);

  %% For each variable find a point.
  points = randi(l - 1, N, var_count);
  flags = toFlag(points);
  
  result = combineWithMask(a, b, flags, l);
end

function h = multiPoint(n)
	 %MULTIPOINT Return a function that produces MULTIPOINT_(N, A, B,
	 % L) when given A, B and L.
	 % If N is 1, it returns SINGLEPOINT.
	 %   H = MULTIPOINT(N)
	 %
	 % See also CROSSOVER>SINGLEPOINT, CROSSOVER>MULTIPOINT_.
  
  if (n == 1)
	h = @singlePoint;
  else
	h = @(a, b, l) multiPoint_(n, a, b, l);
  end
end

function result = multiPoint_(n, a, b, l)
   %MULTIPOINT_ Choose N random points in [1, L - 1] and split both A
   % and B into N+1 parts, then swap and merge them to create two
   % children.
   % N different points are choosen for each variable in A and B.
   %   CHILDREN = MULTIPOINT_(N, A, B, L) Create two children from A
   % and B, with a chromosome length of L.
   %
   % See also CROSSOVER>SINGLEPOINT, CROSSOVER>MULTIPOINT.
  
  global UTILS;
  
  if (n >= l)
	error('multipoint crossover: N must be < L');
  end
  
  dim = size(a);
  N = dim(1);
  var_count = dim(2);
  
  %% NOTE: According to the slides, it should be left to the user to
  %% specify weither or not variables are combined (and l may be
  %% different for each).
  %% This is _not_ handled (btw).
  %% We need the points to be sorted. But they are also random values!
  %% The fastest way to do that (other than using a for-loop) is to
  %% use the indices of those values instead of the values themselves.
  %% We first generate as many random values as the maximum point we
  %% need, then sort the array and use the sorted indices as random
  %% values.
  %% @perf: We generate l - 1 points even if we only need 2...
  %% TODO(@perf): See if using randperm is faster.
  BY_COLUMN = 2;
  [~, indices] = sort(rand(N, l - 1, var_count), BY_COLUMN);
  points = sort(indices(:, 1:n, :), BY_COLUMN);

  
  flags = toFlag(points);

  %% Each point has a corresponding flag. To make the computation
  %% faster, we first compute the final flag (the mask) we get after
  %% splitting at each different point.
  %% This is done by xoring each flag successively.
  %% e,g: points = [2, 4] => [0..00011, 0..01111] => 0..01100.
  mask = zeros(N, var_count);
  for i = 1:var_count
    mask(:, i) = UTILS.reduce(@bitxor, flags(:, :, i), 0);
  end
  
  result = combineWithMask(a, b, mask, l);
end

%% NOTE/FIXME: If t(a) and t(b) is negative and p is a relative sum,
%% the result is bogus! Do not do that!
function h = uniform(p, t)
   %UNIFORM Return a function that produces UNIFORM_(P, A, B, L) when
   % given A, B and L.
   % P Can either be:
   % - A number in [0, 1] (T must not be defined)
   %   If P == 0.5, this calls an optimized version of UNIFORM_:
   %   UNIFORM_05_.
   % - A function handle that takes two parameters and produces a
   %   number in [0, 1].
   %   In that case, T must also be a function handle and, given an
   %   individual, must return the associated value that is then given
   %   to P.
   %   i.e, When given A, B and L, this produces
   %        UNIFORM_(P(T(A), T(B)), A, B, L).
   %
   % See also CROSSOVER>UNIFORM_, CROSSOVER>UNIFORM_05_.
  
  global UTILS;
  
  if (UTILS.isMatlab)
	is_function_handle = @(f) isa(f, 'function_handle');
  end
  
  if (isnumeric(p))
	if ((p < 0) || (p > 1))
	  error('uniform: P must be in [0, 1]');
	elseif (p == 0.5)
	  h = @uniform_05_;
	else
	  h = @(a, b, l) uniform_(p, a, b, l);
	end
  elseif (is_function_handle(p) && is_function_handle(t))
	h = @(a, b, l) uniform_(p(t(a), t(b)), a, b, l);
  else
	error('uniform: either P must be a real in [0, 1] or P and T have to be two function handles');
  end
end

function result = uniform_05_(a, b, l)
   %UNIFORM_05_ This is an optimized version of UNIFORM_(P, A, B, L),
   % when P == 0.5.
   %
   % See also CROSSOVER>UNIFORM, CROSSOVER>UNIFORM_.
  
  %% NOTE: According to the slides, it should be left to the user to
  %% specify weither or not variables are combined (and l may be
  %% different for each).
  %% This is _not_ handled (btw).
  dim = size(a);
  N = dim(1);
  var_count = dim(2);
  max_val = 2^l - 1;

  %% When P == 0.5, every allele has a 50% chance to crossover. This
  %% means that, for a given length l, each crossover pattern (a mask)
  %% has a (1/2)^l chance of happening.
  %% Which is the same as the chance of obtaining a given random
  %% integer in [0, 2^l - 1].
  %% So we can just represent the mask as a random number in this same
  %% interval.
  %% NOTE(@perf): This does not depend on l!
  mask = randi(max_val, N, var_count);
  result = combineWithMask(a, b, mask, l);
end

function result = uniform_(p, a, b, l)
	%UNIFORM_ Create two children from A and B.
	% The first child has a possibility P that a given allele comes
	% from A, and (1 - P) that it comes from B.
	% Those probabilities are inverted for the second child.
	%   CHILDREN = UNIFORM_(P, A, B, L)
	%
	% See also CROSSOVER>UNIFORM, CROSSOVER>UNIFORM_05_.
  
  global UTILS;

  if (any(p < 0 | p > 1))
	error('uniform: P must be in [0, 1]');
  end

  %% NOTE: According to the slides, it should be left to the user to
  %% specify weither or not variables are combined (and l may be
  %% different for each).
  %% This is _not_ handled (btw).
  dim = size(a);
  N = dim(1);
  var_count = dim(2);

  %% For each allele (and for each variable), we generate a random
  %% number in [0, 1] and compare it to p.
  %% This gives us an array of 1s and 0s that we then convert to an
  %% integer (the mask).
  %% For the first child, 1s correspond to the alleles that come from
  %% A and 0s to the ones that come from B.
  %% This is inverted for the second child.
  %%
  %% 'reshape' is used to handle different masks for each variable (we
  %% have an array with Nxlxvar_count dimensions, with N the number of
  %% pairs, l the chromosome length and var_count the number of
  %% variables).
  mask_as_array = rand(N, l, var_count) <= p;
  mask = reshape(UTILS.arrayToDec(mask_as_array), [], var_count);
  
  result = combineWithMask(a, b, mask, l);
end

function result = toFlag(point)
	%TOFLAG Convert POINT into a decimal flag: everything before that
	% point (right-wise) is 1, and everything after is 0.
	% e,g: point = 3 => 0..00111
  
  result = (2 .^ point) - 1;
end

function result = combineWithMask(a, b, mask, l)
	   %COMBINEWITHMASK Create two children based on MASK (an integer
	   % between 0 and 2^L -1).
	   % For the first child, 1s indicate which parts of parent A are
	   % kept. Same goes for 0s and B.
	   % 1s and 0s are inverted for the second child.
  
  max_val = 2^l - 1;
  inv_mask = bitxor(mask, max_val);
  
  result = [ bitor(bitand(a, mask), bitand(b, inv_mask)), ...
			 bitor(bitand(b, mask), bitand(a, inv_mask)) ];
end

%% Arithmetic crossovers
function result = wholeArithmetic(a, b, ~)
	%WHOLEARITHMETIC Create two children from A and B using a linear
	% interpolation (1 - t) * p1 + p2 * t with the same t (a random
	% value in [0, 1]) for each variable (with p1 being either A or B,
	% and p2 being the other).
	%
	% See also CROSSOVER>LOCALARITHMETIC, CROSSOVER>ARITHMETICCROSSOVER_.
  
  result = arithmeticCrossover_(1, a, b);
end

function result = localArithmetic(a, b, ~)
   %LOCALARITHMETIC Create two children from A and B using a linear
   % interpolation (1 - t) * p1 + p2 * t with a different t (a random
   % value in [0, 1]), for each variable (with p1 being either A or B,
   % and p2 being the other).
   %
   % See also CROSSOVER>WHOLEARITHMETIC, CROSSOVER>ARITHMETICCROSSOVER_.
  
  dim = size(a);
  var_count = dim(2);

  result = arithmeticCrossover_(var_count, a, b);
end

function result = arithmeticCrossover_(n, a, b)
	 %ARITHMETICCROSSOVER_ Create two children from A and B using a
	 % linear interpolation (1 - t) * p1 + p2 * t with either the same
	 % or a different t (a random value in [0, 1]), for each variable
	 % (with p1 being either A or B, and p2 being the other).
	 %
	 % See also CROSSOVER>WHOLEARITHMETIC, CROSSOVER>LOCALARITHMETIC.

  %% NOTE/IMPORTANT: n must either be 1 or equal to var_count. Because
  %% this function is not exported, this check is never done, but keep
  %% that in mind ;)
  dim = size(a);
  alpha = rand(dim(1), n);
  beta = 1 - alpha;

  %% Linear interpolation
  result = [(a .* alpha) + (beta .* b), (b .* alpha) + (beta .* a)];
end

function h = blend(alpha)
%BLEND Return a function that produces BLEND_(ALPHA, A, B, CONTEXT)
% when given A, B and CONTEXT.
% ALPHA defaults to 0.5.
%   H = BLEND(ALPHA)
%
% See also CROSSOVER>BLEND_.
  if ~exist('alpha', 'var')
	alpha = 0.5;
  end
  
  h = @(a, b, cx) blend_(alpha, a, b, cx);
end

function result = blend_(alpha, a, b, context)
%BLEND_ Create two children from A and B using ALPHA as the range
% expansion parameter.
%
% Given A(i) and B(i) and assuming A(i) < B(i), produces
% C(i) = rand((A(i) - ALPHA * (B(i) - A(i))), (B(i) + ALPHA * (B(i) - A(i))))
%
% This uses a clamping function.
%
% See also CROSSOVER>BLEND, CLAMP.
  
  constraints = context.constraints;
  
  dim = size(a);
  N = dim(1);
  var_count = dim(2);

  %% The formula assumes a(i) < b(i). But this may not be the case!
  %% For each variable (and each pair), we get the min and max values.
  %% (So we can just replace a(i) by min_vals(i), and b(i) by
  %% max_vals(i) in the formula).
  max_vals = max(a, b);
  min_vals = min(a, b);

  %% The 'ALPHA * (B(i) - A(i))' part
  delta = alpha * (max_vals - min_vals);

  %% Lower and upper bounds of rand
  lb = min_vals - delta;
  ub = max_vals + delta;

  lowest = constraints(:, 1)';
  biggest = constraints(:, 2)';

  clamp_fn = context.clamp_fn;

  result = [ blendChild(lowest, biggest, lb, ub, N, var_count, clamp_fn), ...
			 blendChild(lowest, biggest, lb, ub, N, var_count, clamp_fn) ];
end

function result = blendChild(lowest, biggest, lb, ub, N, var_count, clamp_fn)
  %BLENDCHILD Return N children, with VAR_COUNT variables in [lb, ub].
  % CLAMP_FN is used to clamp out of bound variables inside
  % [LOWEST, BIGGEST].
  
  result = (ub - lb) .* rand(N, var_count) + lb;
  result = clamp_fn(result, lowest, biggest);
end

function h = simulatedBinary(n)
%SIMULATEDBINARY Return a function that produces SIMULATEDBINARY_(N, A, B)
% when given A and B.
% N must be >= 0.
%   H = BLEND(N)
%
% See also CROSSOVER>SIMULATEDBINARY_.
  if (n < 0)
	error('simulatedBinary: N must be >= 0');
  end
  
  h = @(a, b, cx) simulatedBinary_(n, a, b, cx);
end

function result = simulatedBinary_(n, a, b, ~)
%SIMULATEDBINARY_ Create two children from A and B using N to generate
% a random number.
%
% See also CROSSOVER>SIMULATEDBINARY.
  
  dim = size(a);
  N = dim(1);
  var_count = dim(2);

  %% Generate u for all variables (and each pair)
  u = rand(N, var_count);

  %% Sharp S is defined differently if u > 0.5 or u <= 0.5.
  %% To avoid adding a conditional (and therefore a loop), we compute
  %% a logical vector that is 1 if u <= 0.5 and 0 otherwhise.
  %% We then proceed to compute the two different methods using u as a whole.
  %% We use the logical vector to know which parts of those methods
  %% belong to which u (by multiplying by the vector, if u <= 0.5, we
  %% keep the value, otherwhise, we discard it. The opposite is true
  %% if we multiply by the inverse vector).
  below = u <= 0.5;
  
  below_sharp_s = (2 * u) .^ (1 / (n + 1)) .* below;
  above_sharp_s = (2 * (1 - u)) .^ (- (1 / n) + 1) .* ~below;
  
  sharp_s = below_sharp_s + above_sharp_s;
  
  common_part = 0.5 * (a + b);
  delta = 0.5 * sharp_s .* (a - b);

  result = [common_part + delta, common_part - delta];
end

function h = oneBitAdaptation(fn_on_0, fn_on_1)
 %ONEBITADAPTATION Return a function that produces
 % oneBitAdaptation_(FN_ON_0, FN_ON_1, A, B, L) when given A, B and L.
 %
 % See also CROSSOVER>ONEBITADAPTATION_.
  
  h = @(a, b, l) oneBitAdaptation_(fn_on_0, fn_on_1, a, b, l);
end

function result = oneBitAdaptation_(fn_on_0, fn_on_1, a, b, l)
  %ONEBITADAPTATION_ Compare the last bit of A and B. If they are the
  % same, use the associated crossover function (FN_ON_0 if 0, FN_ON_1
  % if 1).
% If they are not the same, a random number in [0, 1] determines which
% one to use.
%
% See also CROSSOVER>ONEBITADAPTATION.
  
  global UTILS;
  
  [N, var_count] = size(a);

  %% Bit xor all variables in a and b so the check of the last bit
  %% depends on every variable (not just the last one).
  %% (NOTE: bitxor is used because it is the only 2-operands bit
  %% function that produces 1 and 0 half the time)
  bit_xor_all = @(x) UTILS.reduce(@bitxor, x, 0);

  last_bit_a = bitand(bit_xor_all(a), 1);
  last_bit_b = bitand(bit_xor_all(b), 1);

  last_bit_diff = logical(bitxor(last_bit_a, last_bit_b));
  last_bit_same = ~last_bit_diff;
  
  result = zeros(N, var_count * 2);

  evalSwitch = generateChoiceEval(fn_on_0, fn_on_1, a, b, l);
  
  last_bit_is_0 = ~last_bit_a(last_bit_same);
  result(last_bit_same, :) = evalSwitch(last_bit_same, last_bit_is_0);
  
  random_bit_value = round(rand(nnz(last_bit_diff), 1));
  result(last_bit_diff, :) = evalSwitch(last_bit_diff, random_bit_value);
end

function h = generateChoiceEval(fn_on_0, fn_on_1, a, b, l)
    h = @(indices, choice) choiceEval_(fn_on_0, fn_on_1, indices, choice, a, b, l);
end

function result = choiceEval_(fn_on_0, fn_on_1, indices, choice, a, b, l)
    result = fn_on_0(a(indices, :), b(indices, :), l) .* choice + ...
             fn_on_1(a(indices, :), b(indices, :), l) .* ~choice;
end
