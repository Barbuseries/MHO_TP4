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


PROFILING = 0;

if (PROFILING)
  profile off;
  profile clear;
  profile on;
end

algo = GA.create(NSGA2);


p = PROBLEM.kursawe(algo, 3);

config = algo.defaultConfig();
config.Pc = 0.9;
config.Pm = 1/52;
config.N = 100;
config.M = 40;
config.G_max = 250;
config.l = 52;
config.crossover_fn = CROSSOVER.singlePoint;
config.mutation_fn = MUTATION.bitFlip;
%config.stop_criteria_fn = STOP_CRITERIA.meanChangeRate(0.005);

r = p.optimize(config);
disp(r);

GA.showPaleto(p, r);

if (PROFILING)
  if (UTILS.isMatlab)
	profile viewer;
  else
	profshow;
  end
end
