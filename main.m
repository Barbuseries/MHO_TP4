%% Includes
Utils;         global UTILS;

Selection;     global SELECTION;
Crossover;     global CROSSOVER;
Mutation;      global MUTATION;
StopCriteria;  global STOP_CRITERIA;
Clamp;         global CLAMP;
Spea2;         global SPEA2;
Problem;       global PROBLEM;


PROFILING = 0;

if (PROFILING)
  profile off;
  profile clear;
  profile on;
end

p = PROBLEM.kursawe(SPEA2, 2);

config = SPEA2.defaultConfig();
config.Pc = 0.8;
config.Pm = 0.01;
config.N = 100;
config.M = 150;
config.G_max = 200;
config.l = 52;
config.crossover_fn = CROSSOVER.uniform(0.5);
config.mutation_fn = MUTATION.bitFlip;

r = p.optimize(config);
disp(r);

SPEA2.showPaleto(p, r);

if (PROFILING)
  if (UTILS.isMatlab)
	profile viewer;
  else
	profshow;
  end
end
