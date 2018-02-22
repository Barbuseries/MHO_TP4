function Selection
	%SELECTION All selection functions.
	%
	% wheel
	% stochasticUniversalSampling
	% tournament(K), K in [1, N]
	% unbiasedTournament(K), K in [1, N]
	% truncation(C), C in [1, N]
	%
	% See also SELECTION>WHEEL, SELECTION>STOCHASTICUNIVERSALSAMPLING,
	% SELECTION>TOURNAMENT, SELECTION>UNBIASEDTOURNAMENT,
	% SELECTION>TRUNCATION
  
  global SELECTION;
  
  SELECTION.wheel = @wheel;
  SELECTION.stochasticUniversalSampling = @stochasticUniversalSampling;
  SELECTION.tournament = @tournament;
  SELECTION.unbiasedTournament = @unbiasedTournament;
  SELECTION.truncation = @truncation;
end

function result = wheel(probabilities)
 %WHEEL For as many times as there are elements in PROBABILITIES, find
 % the first index for which cumsum(PROBABILITIES) >= rand.
 %
 % See also SELECTION>STOCHASTICUNIVERSALSAMPLING.
  
  global UTILS;

  %% We need to select as many individuals as there already are.
  wheel = rand(length(probabilities), 1);
  
  result = UTILS.select(probabilities, wheel);
end

function result = stochasticUniversalSampling(probabilities)
%STOCHASTICUNIVERSALSAMPLING Same as a wheel selection, but instead of
% generating N random numbers, use N equidistant pointers which are
% used as thresholds (the first pointer position is a random number in
% [0, 1/N]), with N being the number of elements in PROBABILITIES.
%
% See also SELECTION>WHEEL.
  
  global UTILS;
  
  N = length(probabilities);

  delta = 1/N; %% Distance between two pointers
  start_pos = delta * rand();
  pointers = delta * (0:(N-1)) + start_pos; %% Generate equidistant pointers starting at start_pos

  %% Find the first index for which
  %% cumsum(probabilities)(i) >= pointers(i)
  result = UTILS.select(probabilities, pointers');
end

function h = tournament(k, c)
	   %TOURNAMENT Return a function that produces TOURNAMENT_(K,
	   % PROBABILITIES, C) when given PROBABILITIES.
	   %   H = TOURNAMENT(K, C = N)
	   %
	   % 1 <= K <= N, with N = length(PROBABILITIES)
	   %
	   % See also SELECTION>UNBIASEDTOURNAMENT, SELECTION>TOURNAMENT_.
  
  if (k < 1)
	error('K must be in [1, N]');
  end
  
  if ~exist('c', 'var')
      h = @(p) tournament_(k, p, length(p));
  else
      h = @(p) tournament_(k, p, c);
  end
end

function result = tournament_(k, probabilities, count)
%TOURNAMENT_ Select K elements in PROBABILITIES at random and keep the index of
% the maximum value COUNT times.
%
% See also SELECTION>TOURNAMENT, SELECTION>UNBIASEDTOURNAMENT.
  
  N = length(probabilities);

  if (k > N)
	error('K must be in [1, N]');
  end

  %% For each selection, select  which k elements we compare.
  random_indices = randi(N, count, k);

  %% For each selection, the index of the maximum value (in random_indices)
  max_indices = tournamentSelect_(random_indices, probabilities);
  result = random_indices(max_indices);
end

function h = unbiasedTournament(k)
%UNBIASEDTOURNAMENT Return a function that produces UNBIASEDTOURNAMENT_(K,
% PROBABILITIES) when given PROBABILITIES.
%   H = UNBIASEDTOURNAMENT(K)
%
% 1 <= K <= N, with N = length(PROBABILITIES)
%
% See also SELECTION>TOURNAMENT, SELECTION>UNBIASEDTOURNAMENT_.
  
  if (k < 1)
	error('K must be in [1, N]');
  end
  
  h = @(p) unbiasedTournament_(k, p);
end

function result = unbiasedTournament_(k, probabilities)
 %UNBIASEDTOURNAMENT_ Create K permutations of probabilities. At each
 % index, compare the K permutations and keep the index of the maximum
 % associated element in PROBABILITIES.
 %
 % See also SELECTION>UNBIASEDTOURNAMENT, SELECTION>TOURNAMENT.
  
  N = length(probabilities);

  if (k > N)
	error('K must be in [1, N]');
  end

  permutations = zeros(k, N);
  for i = 1:k
	permutations(i, :) = randperm(N);
  end

  %% Comparing the permutations column-wise is the same as transposing
  %% them and comparing them row-wise (which is the same as a
  %% tournament selection).
  permutations = permutations';
  max_indices = tournamentSelect_(permutations, probabilities);
  
  result = permutations(max_indices);
end

function result = tournamentSelect_(random_indices, probabilities)
  [N, k] = size(random_indices);
  
  %% If we are only selecting one individual, we do not need to compare
  %% the probabilities.
  if (k == 1)
      result = 1:N;
      return;
  end
  
  %% Find the indices of the maximum values row-wise.
  BY_ROW = 2;
  [~, rel_max_indices] = max(probabilities(random_indices), [], BY_ROW);
  
  %% max returns the index relative to the row (not the whole matrix)
  %% i.e, 1 correponds to the first column, no matter which row we are
  %% in.
  %% And I want 1 to only refer to the first column in the first row.
  result = relativeToExactIndex_(rel_max_indices, N);
end

function result = relativeToExactIndex_(ind, N)
  
  %% ind is the relative index of the column (in 1:N).
  %% What we want is its exact location in the matrix.
  %% We know the matrix has N rows and that matrix indexing is
  %% column-wise.
  %% i.e, 1 refers to (1, 1), 2 to (1, 2), ..., N + 1 to (2, 1), ...
  result = (ind - 1) * N + (1:N)';
end

function h = truncation(c)
			%TRUNCATION Return a function that produces TRUNCATION_(C,
			% PROBABILITES) when given PROBABILITIES.
			%
			% See also SELECTION>TRUNCATION_.
  
  if (c <= 0)
	error('C must be > 0');
  end

  if (c == 1) %% Keep the whole population.
	%% Small optimization, because I can (even though this is not a
	%% real use case...).
	h = @(p) 1:length(p);
  else
	h = @(p) truncation_(c, p);
  end
end

function result = truncation_(c, probabilities)
 %TRUNCATION_ Select N elements by keeping C times the N/C best elements.
 % i.e, C = 10: keep the tenth of the population.
 %
 % If R = mod(N, C) ~= 0, fill the R remaining elements with the first
 % R best elements.
 %
 % See also SELECTION>TRUNCATION.
  
  N = length(probabilities);
  
  [~, ordered_indices] = sort(probabilities, 'descend');

  keep_count = floor(N/c);
  kept_indices = ordered_indices(1:keep_count);

  %% We use the same elements c times.
  fill_count = c * keep_count;

  result = zeros(1, N);
  result(1:fill_count) = repmat(kept_indices, 1, c);

  rest = N - fill_count;

  %% N is not divisible by C.
  %% Get as many elements from kept_indices as needed (rest <
  %% length(kept_indices)).
  if (rest ~= 0)
	result(fill_count+1:end) = kept_indices(1:rest);
  end
end
