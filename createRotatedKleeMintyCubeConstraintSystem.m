function [problem] = createRotatedKleeMintyCubeConstraintSystem(dimension)
  % [A, b] = createKleeMintyCubeConstraintSystem(dimension,angle)
  %    Creates a linear system of equations
  %        Ax <= b
  %    that represents the feasible region.
  %
  %    The constraints represent the Klee-Minty cube 
  %    Klee, V., Minty, G.J.:
  %    How good is the simplex algorithm? In: Shisha, O. (ed.) Inequalities III,
  %    pp. 159–175. Academic, New York (1972)
  %
  %    in the representation from 
  %
  %    Deza, A., Nematollahi, E., Peyghami, R., Terlaky, T.: 
  %    The central path visits all the vertices of the klee–minty cube. 
  %    Optimisation Methods and Software 21(5), 851--865 (2006)
  %
  %    The actual constraints read
  %    
  %                    0 <= x_{1} <= 1
  %    EPSILON * x_{k-1} <= x_{k} <= 1 - EPSILON * x_{k-1}  forall k = 2,...,dimension 
  %
  %    with positive real number EPSILON < 1/3. 
  %
  %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %    The objective is to minimize the linear function 
  %    f(x) = c' * x = x_{dimension}
  %    the coefficients of f are returned in the vector
  %             c = (0, 0, 0, ..., 0, 1)
  %
  %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %    OUTPUT: 
  %         problem.A       - matrix of the inequality constraints
  %         problem.b       - right side of the linear constraints 
  %         problem.c       - vector of the linear objective function coefficients
  %         problem.R       - rotation matrix of specified angle
  %         problem.t       - translation vector
  %
  %         problem.Fname   - objective function name
  %
  %         problem.lower_bounds  - problem specific lower parameter vector bounds
  %         problem.upper_bounds  - upper parameter vector bounds
  %
  %         as well as information on the optimal solutions and target
  %         definitions for postprocessing purposes only (see below)
  
  
  A = [];
  b = [];
  c = [];

  if (dimension <= 0)
    error('dimension must be greater than zero');
  end

  epsi = 1/10;
  
  A1 = zeros(dimension, dimension);
  A2 = zeros(dimension, dimension);
  b = ones(dimension, 1);
  c = zeros(dimension, 1);
  A1(1,1) = 1;
  A2(1,1) = -1;
  for i = 2:dimension
    A1(i, i)  = 1;
    A1(i,i-1) = epsi;
    A2(i, i)  = -1;
    A2(i,i-1) = epsi;
  end
  
  c(dimension,1) = 1;
  
  problem.dim = dimension;
  
  % incorporation of the component-wise non-negativity constraints  
  problem.A     = [A1;A2];
  problem.b     = [b;zeros(dimension,1)];
  problem.c     = c;
  
  % Generate motion (rotation + translation) for rotated Klee-Minty problem
  angle         = 10*pi/180;           % rotation angle in radian measure
  
  rotv1         = [zeros(dimension-1,1);1];
  rotv1         = rotv1./norm(rotv1);   % normalized rotation axis 1
  rotv2         = [ones(dimension-1,1);0];
  rotv2         = rotv2./norm(rotv2);   % normalized rotation axis 2
  
  % The rotation is performed in the hyperplane defined by rotv1 and rotv2
  V             = rotv1*rotv1' + rotv2*rotv2';
  W             = rotv1*rotv2' - rotv2*rotv1';
  E             = eye(dimension);
  
  problem.R     = (E  + (cos(angle)-1).*V - sin(angle).*W);
  problem.t     = ones(dimension,1).*dimension^3;
  
  % box-constraints (for optional use)
  problem.lower_bounds = -ones(dimension,1).*0;
  problem.upper_bounds = ones(dimension,1).*5*dimension^3;
  
  % objective function
  problem.Fname = 'evaluateRotatedKleeMintyProblem';
  
  % optimal parameter vector
  problem.yOpt  = [zeros(1,dimension)]';
  problem.xOpt  = problem.R'*problem.yOpt+problem.t;
  problem.consumed = 0;
  
  % optimal objective function value
  eval(['problem.fOpt =' problem.Fname '(problem.xOpt,problem);']);
  
  % create targtes for runtime measurements
  
    problem.ConTarNum   = 51;
    problem.ConGp       = 4;
    problem.ConLp       = -6;
    problem.ConTarInt   = [problem.ConGp:(problem.ConLp-problem.ConGp)/(problem.ConTarNum-1):problem.ConLp];
    problem.ConTarDev   = [10.^problem.ConTarInt , 0];
               
    problem.FunTarNum   = 51;
    problem.FunGp       = 0;
    problem.FunLp       = -8;
    problem.FunTarInt   = [problem.FunGp:(problem.FunLp-problem.FunGp)/(problem.FunTarNum-1):problem.FunLp];
    problem.FunTarDev   = 10.^problem.FunTarInt;
    problem.Ref_Val     = problem.fOpt;
    
    problem.TarNum      = problem.ConTarNum +problem.FunTarNum +1;
  
end
