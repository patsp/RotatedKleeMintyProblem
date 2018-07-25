function [f,g] = evaluateRotatedKleeMintyProblem(x,problem)
    global target_flag
    global consumed
    
    y = problem.R*(x-problem.t);
    f = sum(problem.c.*x); %% y
    g = problem.A*y-problem.b;
    
    if target_flag > 0
        consumed = consumed+1;
        logRotatedKleeMintyEvals(problem,f,g);
    end
    
end
