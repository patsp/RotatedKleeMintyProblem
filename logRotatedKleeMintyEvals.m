function [] = logRotatedKleeMintyEvals(problem,f,gv,x)
    global target_flag
    global T
    global Targets
    global consumed
    global best
    global global_best_log
    
    Target_Num      = problem.ConTarNum +problem.FunTarNum + 1; 
    
    if target_flag==1
        T               = [];
        Targets         = zeros(Target_Num,4);
        Targets(:,1:2)  = [[zeros(problem.ConTarNum+1,1),[problem.ConTarDev]'];[(problem.FunTarDev)', zeros(problem.FunTarNum,1)]];
        best            = [f, sum(gv.*(gv>0))];
        target_flag     = 2;
        global_best_log = x;
    end
    
    if ( best(2)==0 && sum(gv.*(gv>0))==0 ) 
        if best(1)>=f
            best = [f, sum(gv.*(gv>0))];
            global_best_log = x;
        end
    elseif  best(2)>= sum(gv.*(gv>0))
        best = [f, sum(gv.*(gv>0))];
        global_best_log = x;
    end
    
    T = [T;[consumed, best(1), best(2)]];
    
    for k=1:Target_Num
       if Targets(k,4)==0
           if k<=problem.ConTarNum+1
                if best(2)<=problem.ConTarDev(k)
                   Targets(k,3) = consumed;
                   Targets(k,4) = 1;
                end
           else
                if abs(problem.Ref_Val-best(1))<=problem.FunTarDev(k-(problem.ConTarNum+1)) && best(2)==0
                   Targets(k,3) = consumed;
                   Targets(k,4) = 1;
                end
           end
       end
    end
end

