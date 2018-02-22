%% IMPORTANT: To make functions work with plots, they must accept
%% matrices as variables (i,e., x = [x11, x12, ...; x21, x22, ...;
%% ...], y = [y11, y12, ...; y21, y22, ...; ...]).
%% Remember that when writing them!

function Problem
  global PROBLEM;
  
  PROBLEM.TOTO = @TOTO;
  PROBLEM.kursawe = @kursawe;
  PROBLEM.shaffer = @shaffer;
  PROBLEM.fonsecaFlemming = @fonsecaFlemming;
  PROBLEM.poloni = @poloni;
  PROBLEM.zdt1 = @zdt1;
end

%% TOTO
function result = TOTO(ga)
  result.objective_vector = {@TOTO_f1_, @TOTO_f2_};
  result.constraints = [[0, 10]
					    [0, 10]];

  result.optimize = optimize_(ga, result, 0);
end

function result = TOTO_f1_(x, y)
  result = x - 3 * y;
end

function result = TOTO_f2_(x, y)
  result = 3 * x + y;
end

%% Kursawe
function result = kursawe(ga, n)
  result.objective_vector = {generate_fn_n_(n, @kursawe_f1_), generate_fn_n_(n, @kursawe_f2_)};
  result.constraints = repmat([-5, 5], n, 1);
  
  result.optimize = optimize_(ga, result, 0);
end

function result = kursawe_f1_(n, varargin)
  xi = shapeVariables(n, varargin{:});
  
  xi_sq = xi(:, :, 1:end-1) .^ 2;
  xi_plus_one_sq = xi(:, :, 2:end) .^ 2;
  
  BY_DEPTH = 3;
  result = sum(-10 * exp(-0.2 * sqrt(xi_sq  + xi_plus_one_sq)), BY_DEPTH);
end

function result = kursawe_f2_(n, varargin)
  xi = shapeVariables(n, varargin{:});
  
  BY_DEPTH = 3;
  result = sum(abs(xi).^0.8 + 5*sin(xi.^3), BY_DEPTH);
end


%% Shaffer
function result = shaffer(ga)
  result.objective_vector = {@shaffer_f1_, @shaffer_f2_};
  result.constraints = [-10^3, 10^3];

  result.optimize = optimize_(ga, result, 0);
end

function result = shaffer_f1_(x)
  result = x .^ 2;
end

function result = shaffer_f2_(x)
  result = (x  - 2) .^ 2;
end

%% Fonseca Flemming
function result = fonsecaFlemming(ga, n)
  result.objective_vector = {generate_fn_n_(n, @fonsecaFlemming_f1_), generate_fn_n_(n, @fonsecaFlemming_f2_)};
  result.constraints = repmat([-4, 4], n, 1);

  result.optimize = optimize_(ga, result, 0);
end

function result = fonsecaFlemming_f1_(n, varargin)
  xi = shapeVariables(n, varargin{:});

  BY_DEPTH = 3;
  inv_sqrt3 = 1 / sqrt(3);
  result = 1 - exp(-sum((xi - inv_sqrt3) .^2, BY_DEPTH));
end

function result = fonsecaFlemming_f2_(n, varargin)
  xi = shapeVariables(n, varargin{:});

  BY_DEPTH = 3;
  inv_sqrt3 = 1 / sqrt(3);
  result = 1 - exp(-sum((xi + inv_sqrt3) .^2, BY_DEPTH));
end

%% Poloni
function result = poloni(ga)
  result.objective_vector = {@poloni_f1_, @poloni_f2_};
  result.constraints = [[-pi, pi],
						[-pi, pi]];

  result.optimize = optimize_(ga, result, 0);
end

function result = poloni_f1_(x, y)
  A1 = poloni_b1_(1, 2);
  A2 = poloni_b2_(1, 2);

  result = 1 + ((A1 - poloni_b1_(x, y)) .^2) + ((A2 - poloni_b2_(x, y)) .^2);
end

function result = poloni_f2_(x, y)
  result = ((x + 3) .^2) + ((y + 1) .^ 2);
end

function result = poloni_b1_(x, y)
  result = 0.5 * sin(x) - 2 * cos(x) + sin(y) - 1.5 * cos(y);
end

function result = poloni_b2_(x, y)
  result = 1.5 * sin(x) - cos(x) + 2 * sin(y) - 0.5 * cos(y);
end


%% ZDT1
function result = zdt1(ga, n)
  result.objective_vector = {generate_fn_n_(n, @zdt1_f1_), generate_fn_n_(n, @zdt1_f2_)};
  result.constraints = repmat([0, 1], n, 1);

  result.optimize = optimize_(ga, result, 0);
end

function result = zdt1_f1_(n, varargin)
  xi = shapeVariables(n, varargin{:});

  result = xi(:, :, 1);
end

function result = zdt1_f2_(n, varargin)
  xi = shapeVariables(n, varargin{:});

  x1 = xi(:, :, 1);
  g_x = zdtx_g1_(xi, n);
  
  result = g_x .* (1 - sqrt(x1 ./ g_x));
end

function result = zdtx_g1_(x, n)
  BY_DEPTH = 3;
  
  result = 1 + 9 * sum(x(:, :, 2:end), BY_DEPTH) / (n - 1);
end


function result = optimize_(ga, problem, maximize)
    result = @(config) ga.optimize(maximize, problem.objective_vector, problem.constraints, config);
end

function h = generate_fn_n_(n, fn)
  h = @(varargin) fn(n, varargin{:});
end

%% IMPORTANT: To sum variables (x + y, ...), the sum must be done on
%% the _depth_ (i.e, sum(..., 3)).
%% (x is represented by result(:, :, 1), y by result(:, :, 2), ...)
function result = shapeVariables(n, varargin)
  %% Reshape to have separated variables represented by the depth
  %% (v(:, :, i)) (so it works with matrices)
  [N, ~] = size(varargin{1});
  result  = reshape([varargin{:}], [], N, n);
end
