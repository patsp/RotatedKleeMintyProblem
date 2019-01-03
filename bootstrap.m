function FEperTargetBootstrapped = bootstrap(problem, ...
                                             FEperTarget, ...
                                             samplesize, ...
                                             maxTries, ...
                                             maxRunlength)
  numTargets = size(FEperTarget, 1);
  numRuns = size(FEperTarget, 2) - 2;
  % first two columns are the targets
  FEperTargetBootstrapped = zeros(numTargets, 2 + samplesize);
  FEperTargetBootstrapped(:, 1:2) = FEperTarget(:, 1:2);
  % for all targets
  for targetIdx = 1:numTargets
    % if there is at least one success for the current target
    if any(FEperTarget(targetIdx, 3:end) ~= 0)
      for s = 1:samplesize
        % sample till the first successful run is sampled
        runlengthSum = 0;
        success = false;
        numTries = 0;
        while ~success && numTries < maxTries
          randomIdx = randi([3, 3 + numRuns - 1], 1, 1);
          runlength = FEperTarget(targetIdx, randomIdx);
          if runlength == 0
            runlengthSum = runlengthSum + maxRunlength;
          else
            success = true;
            runlengthSum = runlengthSum + runlength;
          end
          numTries = numTries + 1;
        end
        FEperTargetBootstrapped(targetIdx, 3 + s - 1) = runlengthSum;
        % ensure that more difficult target is not reached with
        % fewer evaluations than less difficult target
        if targetIdx > 1
          FEperTargetBootstrapped(targetIdx, 3 + s - 1) = ...
          max(FEperTargetBootstrapped(targetIdx - 1, 3 + s - 1), ...
              FEperTargetBootstrapped(targetIdx, 3 + s - 1));
        end
      end
    else
      % no success possible in the bootstrapped runs
      FEperTargetBootstrapped(targetIdx, 3:end) = 0;
    end
  end
end
