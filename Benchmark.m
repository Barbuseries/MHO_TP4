%% FIXME: Depracated!
%% TODO: Update!
function [v, v2, l] = Benchmark(p, config, count, leg, v, v2, l)
  benchmark(1:count) = struct('h', 0);

  for i = 1:count
    [~, h] = p.optimize(config);
    benchmark(i).h = h;
    disp(i)
  end
  
  a = [benchmark.h];
  b = [a.very_best];
  iteration_counts = cellfun(@length, {a.iterations});
  v(:, end+1) = iteration_counts';
  v2(:, end+1) = [b.fitness]';
  l(end+1) = {leg};
end
