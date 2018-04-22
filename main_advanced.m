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
Pesa;          global PESA;
Pesa2;         global PESA2;
Ibea;          global IBEA;
Problem;       global PROBLEM;
Benchmark;     global BENCHMARK;


PROFILING = 0;

if (PROFILING)
  profile off;
  profile clear;
  profile on;
end

BENCHMARKING = 0;
COMPARING = 1;
RUNNING = 2;


%% When not bencharmking, those can be uncommented (as well as their
%% var count below);

all_problem_handles = {PROBLEM.schaffer, PROBLEM.kursawe, ...
                       PROBLEM.fonsecaFlemming, PROBLEM.poloni, ...
                       PROBLEM.zdt1, PROBLEM.zdt2, ...
                       PROBLEM.zdt3, PROBLEM.zdt4, ...
                       PROBLEM.zdt6};
%% Minus sign: already preconfigured
all_var_counts = [-1, 3, ...
                   3, -2, ...
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

all_algos = [NSGA2, SPEA2, PESA, PESA2, IBEA, IBEA, IBEA];


%% Set this to any of the above.
doing = COMPARING;
if (doing == BENCHMARKING)
    all_configs = cell(length(all_algos), length(all_problems));

    N = 100;
    G_MAX = 1000;
    l = -1;
    CROSSOVER_FN = CROSSOVER.simulatedBinary(20);
    MUTATION_FN = MUTATION.polynomial(20);
    %% CROSSOVER_FN = CROSSOVER.uniform(0.5);
    %% MUTATION_FN = MUTATION.bitFlip;

    %% Common config
    for j = 1:length(all_algos)
        for i = 1:length(all_problems)
            if (l == -1)
              Pm = 1 / abs(all_var_counts(i));
            else
              Pm = 1 / l;
            end

            all_configs{j, i} = struct('N', N, 'l', l, 'Pc', 0.8, 'Pm', Pm,'G_max', G_MAX, ...
									   'crossover_fn', CROSSOVER_FN, 'mutation_fn', MUTATION_FN);
        end
    end

    %% NSGA2 config
    ALGO_INDEX = 1;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% SPEA2 config
    ALGO_INDEX = 2;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.M = 100;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% PESA config
    ALGO_INDEX = 3;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.M = 100;
        config.C = 32;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% PESA2 config (same as PESA)
    ALGO_INDEX = 4;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.M = 100;
        config.C = 32;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% IBEA config
    ALGO_INDEX = 5;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.Pc = 1;
        config.Pm = 0.1;
        config.kappa = 0.002;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% IBEA adaptative (0.002) config
    ALGO_INDEX = 6;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.Pc = 1;
        config.Pm = 0.1;
        config.kappa = 0.002;
        config.adaptative = 1;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    %% IBEA adaptative (0.05) config
    ALGO_INDEX = 7;
    for i = 1:length(all_problems)
        config = all_configs{ALGO_INDEX, i};
        config.Pc = 1;
        config.Pm = 0.1;
        config.kappa = 0.05;
        config.adaptative = 1;

        all_configs{ALGO_INDEX, i} = config;
    end
    
    data = BENCHMARK.run(all_problems, all_algos, all_configs, 30);
    %BENCHMARK.plot(all_algos, all_problems, data, true)
elseif (doing == COMPARING)
    algo_count = length(all_algos);
    problem_count = length(all_problems);
    
    for j = 1:problem_count
        figure(j)
        clf;
        hold on;
        
        fprintf(1, '%d/%d\n', j, problem_count);
        
        all_h = [];
        plot_legend = cell(0);
        
        for i = [1, 2, 3, 6]%%1:algo_count
            algo = GA.create(all_algos(i));

            p = all_problems{j}(algo);

            config = algo.defaultConfig();
            config.l = -1;
            config.Pc = 0.8;
            %%config.Pm = 1 / config.l;
            config.Pm = 1 / abs(all_var_counts(j));
            config.N = 100;
            config.M = 100;
            config.crossover_fn = CROSSOVER.simulatedBinary(20);
            config.mutation_fn = MUTATION.polynomial(20);
            
            if ((i == 3) || (i == 4))
                config.C = 32;
            end
            
            if ((i == 5) || (i == 6) || (i == 7))
               config.Pc = 1;
               config.Pm = 0.1;
               
               if ((i == 5) || (i == 6))
                   config.kappa = 0.002;
               end
               
               if (i == 6)
                   config.adaptative = 1;
                   algo.name = "IBEA+(0.002)";
               end
               
               if (i == 7)
                   config.kappa = 0.05;
                   config.adaptative = 1;
                   algo.name = "IBEA+(0.05)";
               end
            end
            
            config.G_max = 200;
            %config.stop_criteria_fn = STOP_CRITERIA.meanChangeRate(0.005);

            [r, h] = p.optimize(config);
            
            show_optimal_front = ~isempty(p.optimal_solutions(1));
            [all_h, plot_legend] = GA.iterativeShowPaleto(p, r, show_optimal_front, algo.name, all_h, plot_legend);
        end
    end
elseif (doing == RUNNING)
    algo = GA.create(SPEA2);

    p = PROBLEM.schaffer(algo);

    config = algo.defaultConfig();
    config.Pc = 0.9;
    config.Pm = 1 / 5;
    config.N = 100;
    config.M = 100;
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
