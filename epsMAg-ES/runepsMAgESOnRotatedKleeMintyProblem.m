function [out, global_best] = runepsMAgESOnRotatedKleeMintyProblem(problemIn, ...
                                                                   budget, ...
                                                                   lbnds, ...
                                                                   ubnds, ...
                                                                   inputIn)
  problem = problemIn;
  problem.lower_bounds = lbnds(:);
  problem.upper_bounds = ubnds(:);
  numConstraints = size(gHelper(zeros(size(lbnds(:), 1), 1), problemIn), 1);
  problem.gn = @(~) numConstraints;
  problem.hn = @(~) 0;
  problem.constr_fun_name = @evaluateHelper;
  function [fval, gv, hv] = evaluateHelper(x, ~)
    hv = [];
    [fval, gv] = fgHelper(x, problemIn);
  end
  input = struct();
  D = problemIn.dim;
  input.budget            = budget;
  input.maxIter           = floor(budget / 10);
  input.delta             = 10^-4; % error margin for equality constraints
  input.dim               = D;
  input.lambda            = 4*D;                          % population size
  input.sigma             = 1;                            % initial mutation strength
  input.mu                = floor(input.lambda/3);        % parental ppopulation size
  input.weights = log(input.mu+1/2)-log(1:input.mu)';     % muXone array for weighted recombination
  input.weights = input.weights./sum(input.weights);      % normalize recombination weights array
  input.mueff=1/sum(input.weights.^2);                    % variance-effectiveness of sum w_i x_i
  input.cs = (input.mueff+2) / (D+input.mueff+5);         % t-const for cumulation for sigma control
  input.c1 = 2 / ((D+1.3)^2+input.mueff);                 % learning rate for rank-one update of M
  input.cmu = min(1-input.c1, 2 * (input.mueff-2+1/input.mueff) / ((D+2)^2+input.mueff));     % and for rank-mu update of M
  input.damps = 1 + 2*max(0, sqrt((input.mueff-1)/(D+1))-1) + input.cs;                       % damping for sigma usually close to 1
  problem.fOpt = problemIn.fOpt;
  [out, global_best] = epsMAgES(problem, input, -1);
end

function [f, g] = fgHelper(x, problem)
  [f, g] = evaluateRotatedKleeMintyProblem(x(:), problem);
end

function f = fHelper(x, problem)
  [f, ~] = fgHelper(x, problem);
end

function g = gHelper(x, problem)
  [~, g] = fgHelper(x, problem);
end

