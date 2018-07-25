function [ecdf_data] = assessRotatedKleeMintyPerformance(problem,FEperTarget)
      
    FevInt=[problem.dim:problem.dim:problem.budget+problem.dim];
    for k=1:length(FevInt)
        RR(k) = length(find(FEperTarget(:,3:problem.number_of_runs+2)<FevInt(k) & FEperTarget(:,3:problem.number_of_runs+2)~=0) );
    end

    ecdf_data.IntervalFevals    = FevInt;
    ecdf_data.FevalsPerTarget   = FEperTarget;
    ecdf_data.TargetsPerFevals  = RR;
end
