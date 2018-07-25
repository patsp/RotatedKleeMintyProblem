function [fBest,global_best]=RandS(problem,budget,lower_bounds,upper_bounds,input)
%% Implementation of random search as baseline comparison and for demonstration purposes
%
% Initialization:
    dim                 = input.dim;
  	g                   = 0;
    termination         = 0;
    fOpt                = problem.fOpt;

    % Sample first candidate solution uniformly between lower and upper parameter vector bounds 
    new_solution        = lower_bounds ...
                            + (upper_bounds-lower_bounds).*rand(dim,1); 
    % Evaluation of the constrained function
    [f, gv]             = feval(problem.Fname,new_solution,problem);
    
    % Calculate constraint violation
    conv                = sum(gv.*(gv>0));
    
    % Count function evaluations
    evals               = 1;    
        
    % Initialize best solution found so far
    global_best.y       = new_solution;
	global_best.val     = f;
    global_best.conv    = conv;
    
    while ~termination
        % Randomly sample new candidate solution
        new_solution    = lower_bounds...
                            +(upper_bounds-lower_bounds).*rand(dim,1); 
        
        [f, gv]         = feval(problem.Fname,new_solution,problem);
        conv            = sum(gv.*(gv>0));
        evals           = evals + 1;
        
        % Replace best so far solution if improvement is discovered
        % The selection is based on a lexicographic ordering.
        if (conv==global_best.conv && f <= global_best.val) || conv<global_best.conv
            global_best.y    = new_solution;
            global_best.val  = f;
            global_best.conv = conv;
        end
                        
        % Termination after exceeding the budget of function evaluations or
        % after reaching the final target
        if (abs(global_best.val-fOpt) <= 10^problem.FunLp && global_best.conv==0 ) || evals >= budget 
            termination = 1;
        end        
    end
    fBest               = global_best.val ;   
end
