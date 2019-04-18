function finalTargetHit = isRotatedKleeMintyFinalTargetHit(problem)
  global Targets;

  nTargets = problem.ConTarNum + problem.FunTarNum + 1;
  finalTargetHit = (Targets(nTargets, 4) == 1);
end
