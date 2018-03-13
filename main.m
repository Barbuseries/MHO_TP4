folder = fileparts(which(mfilename)); 
addpath(genpath(folder))

%% Includes
Utils;         global UTILS;

Selection;     global SELECTION;
Crossover;     global CROSSOVER;
Mutation;      global MUTATION;
StopCriteria;  global STOP_CRITERIA;
Clamp;         global CLAMP;
Ga;            global GA;
Spea2;         global SPEA2;
Nsga2;         global NSGA2;
Problem;       global PROBLEM;
Benchmark;     global BENCHMARK;


PROFILING = 0;

if (PROFILING)
  profile off;
  profile clear;
  profile on;
end

BENCHMARKING = 0;
COMPARING_PARETO = 1;
RUNNING = 2;


all_problem_handles = {PROBLEM.schaffer, ...%PROBLEM.kursawe, ...
                       PROBLEM.fonsecaFlemming, ...%PROBLEM.poloni, ...
                       PROBLEM.zdt1, PROBLEM.zdt2, ...
                       PROBLEM.zdt3, PROBLEM.zdt4, ...
                       PROBLEM.zdt6};
all_var_counts = [-1, ...% 3,
                   3, ...%-2,
                   30, 30, 30, 10, 10];

all_problems = cell(1, length(all_var_counts));
for i = 1:length(all_var_counts)
    var_count = all_var_counts(i);

    if (var_count < 0)
        all_problems{i} = all_problem_handles{i};
    else
        all_problems{i} = @(ga) all_problem_handles{i}(ga, var_count);
    end
end

all_algos = [NSGA2, SPEA2];


doing = COMPARING_PARETO;
if (doing == BENCHMARKING)
    all_configs = cell(length(all_algos), length(all_problems));

    N = 100;
    G_MAX = 250;
    l = 52;
    CROSSOVER_FN = CROSSOVER.singlePoint;
    MUTATION_FN = MUTATION.bitFlip;

    %% Common config
    for j = 1:length(all_algos)
        for i = 1:length(all_problems)
            if (l == -1)
                Pm = l / abs(all_var_counts(i));
            else
                Pm = 1 / l;
            end

            all_configs{j, i} = struct('N', N, 'l', l, 'Pc', 0.9, 'Pm', Pm, 'G_max', G_MAX, 'crossover_fn', CROSSOVER_FN, 'mutation_fn', MUTATION_FN);
        end
    end

    %% NSGA2 config
    ALGO_INDEX = 1;

    %% SPEA2 config
    ALGO_INDEX = 2;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.M = ceil(config.N * 1 / 4);

        if (mod(config.M, 2) == 1)
            config.M = config.M + 1;
        end

        all_configs{ALGO_INDEX, i} = config;
    end

    data = BENCHMARK.run(all_problems, all_algos, all_configs, 30);
elseif (doing == COMPARING_PARETO)
    algo_count = length(all_algos);
    problem_count = length(all_problems);
    
    for j = 1:problem_count
        figure(j)
        clf;
        hold on;
        
        all_h = [];
        plot_legend = cell(0);
        
        for i = 1:algo_count
            algo = GA.create(all_algos(i));

            p = all_problems{j}(algo);

            config = algo.defaultConfig();
            config.Pc = 0.9;
            config.Pm = 1 / 2;%all_var_counts(j);
            config.N = 100;
            config.M = 100;
            config.G_max = 250;
            config.l = -1;
            config.crossover_fn = CROSSOVER.simulatedBinary(20);
            config.mutation_fn = MUTATION.polynomial(20);
            %config.stop_criteria_fn = STOP_CRITERIA.meanChangeRate(0.005);

            [r, h] = p.optimize(config);
            disp(r);

            show_optimal_front = (i == 1);
            [all_h, plot_legend] = GA.showPaleto(p, r, show_optimal_front, algo.name, all_h, plot_legend);
        end
    end
elseif (doing == RUNNING)
    algo = GA.create(NSGA2);

    p = PROBLEM.fonsecaFlemming(algo, 2);

    config = algo.defaultConfig();
    config.Pc = 0.9;
    config.Pm = 1 / 2;
    config.N = 100;
    config.M = 26;
    config.G_max = 250;
    config.l = -1;
    config.crossover_fn = CROSSOVER.simulatedBinary(20);
    config.mutation_fn = MUTATION.polynomial(20);
    %config.stop_criteria_fn = STOP_CRITERIA.meanChangeRate(0.005);

    [r, h] = p.optimize(config);
    disp(r);

    GA.showPaleto(p, r, true);
end

if (PROFILING)
  if (UTILS.isMatlab)
	profile viewer;
  else
	profshow;
  end
end
