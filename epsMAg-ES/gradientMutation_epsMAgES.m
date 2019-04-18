function [y]=gradientMutation_epsMAgES(problem,x,gg,hh,funcno)
    eta     = 1e-4;                         % hard coded deviation for determining finite differences               
    n       = length(x);
    Dd      = eye(n);
    dx      = repmat(x,1,n)+eta.*Dd;
    for i=1:n
        [fff, gv, hv] = feval(problem.constr_fun_name,dx(:,i)',funcno);
		dg= [gv;hv];
	    dCx(:,i)  = dg;
    end
    deltaG  = max(0,gg);
    Cx= [gg;hh];

    nabC    = 1/eta.*( dCx - repmat(Cx,1,n));
    delC    = [deltaG;hh];
    inv_nabC= pinv(nabC,1e-12);            % Moore-Penrose inverse of nabC
    deltaX  = -inv_nabC*delC;
    y       = (x+deltaX);
end
