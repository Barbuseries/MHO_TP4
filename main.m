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
Ibea;          global IBEA;
Problem;       global PROBLEM;


PROFILING = 0;

if (PROFILING)
  profile off;
  profile clear;
  profile on;
end

algo = GA.create(IBEA);

p = PROBLEM.zdt3(algo, 30);

config = algo.defaultConfig();
config.l = 52;
config.Pc = 0.7;
config.Pm = 1 / config.l;
config.C = 64;
config.N = 100;
config.M = 100;
config.G_max = 200;
%%config.crossover_fn = CROSSOVER.simulatedBinary(20);
%%config.mutation_fn = MUTATION.polynomial(20);
config.crossover_fn = CROSSOVER.uniform(0.5);
config.mutation_fn = MUTATION.bitFlip;
%%config.stop_criteria_fn = STOP_CRITERIA.meanChangeRate(0.005);

[r, h] = p.optimize(config);
disp(r);

GA.showPaleto(p, r, true);

if (PROFILING)
  if (UTILS.isMatlab)
	profile viewer;
  else
	profshow;
  end
end
