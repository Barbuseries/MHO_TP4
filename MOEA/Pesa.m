function Pesa
  global PESA;

  PESA = Pesa_common("Pesa", @individual_based_selection_);
end

function h = individual_based_selection_(population_grid_index, all_grid_indices, all_squeeze_factors, selection_fn)
  [i, ~] = find(population_grid_index == all_grid_indices');
  squeeze_factor = all_squeeze_factors(i);

  h = @(count) selection_fn(-squeeze_factor, count);
end
