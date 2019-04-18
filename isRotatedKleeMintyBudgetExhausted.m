function budgetExhausted = isRotatedKleeMintyBudgetExhausted(problem)
  global consumed;

  budgetExhausted = (consumed >= problem.budget);
end
