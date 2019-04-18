% Main routine for running the CMA-ES variants on the 2017 CEC Benchmarks
clear all, clc
DIM =  [10,30,50,100];       
    for lll=1:4

        global  initial_flag
        initial_flag = 0;

        % choose the problem dimensionality
        D=DIM(lll)

        % bound constraint definitions for all 18 test functions
        Xmin1=-100*ones(1,D);
        Xmax1=+100*ones(1,D);
        Xmin2=-100*ones(1,D);
        Xmax2=+100*ones(1,D);
        Xmin3=-100*ones(1,D);
        Xmax3=+100*ones(1,D);
        Xmin4=-10*ones(1,D);
        Xmax4=+10*ones(1,D);
        Xmin5=-10*ones(1,D);
        Xmax5=+10*ones(1,D);
        Xmin6=-20*ones(1,D);
        Xmax6=+20*ones(1,D);
        Xmin7=-50*ones(1,D);
        Xmax7=+50*ones(1,D);
        Xmin8=-100*ones(1,D);
        Xmax8=+100*ones(1,D);
        Xmin9=-10*ones(1,D);
        Xmax9=+10*ones(1,D);
        Xmin10=-100*ones(1,D);
        Xmax10=+100*ones(1,D);
        Xmin11=-100*ones(1,D);
        Xmax11=+100*ones(1,D);
        Xmin12=-100*ones(1,D);
        Xmax12=+100*ones(1,D);
        Xmin13=-100*ones(1,D);
        Xmax13=+100*ones(1,D);
        Xmin14=-100*ones(1,D);
        Xmax14=+100*ones(1,D);
        Xmin15=-100*ones(1,D);
        Xmax15=+100*ones(1,D);
        Xmin16=-100*ones(1,D);
        Xmax16=+100*ones(1,D);
        Xmin17=-100*ones(1,D);
        Xmax17=+100*ones(1,D);
        Xmin18=-100*ones(1,D);
        Xmax18=+100*ones(1,D);
        Xmin19=-50*ones(1,D);
        Xmax19=+50*ones(1,D);
        Xmin20=-100*ones(1,D);
        Xmax20=+100*ones(1,D);
        Xmin21=-100*ones(1,D);
        Xmax21=+100*ones(1,D);
        Xmin22=-100*ones(1,D);
        Xmax22=+100*ones(1,D);
        Xmin23=-100*ones(1,D);
        Xmax23=+100*ones(1,D);
        Xmin24=-100*ones(1,D);
        Xmax24=+100*ones(1,D);
        Xmin25=-100*ones(1,D);
        Xmax25=+100*ones(1,D);
        Xmin26=-100*ones(1,D);
        Xmax26=+100*ones(1,D);
        Xmin27=-100*ones(1,D);
        Xmax27=+100*ones(1,D);
        Xmin28=-50*ones(1,D);
        Xmax28=+50*ones(1,D);

        % number of constraint functions per problem

        problem.gn=[1 1 1 2 2 0 0 0 1 0 1 2 3 1 1 1 1 2 2 2 2 3 1 1 1 1 2 2];
        problem.hn=[0 0 1 0 0 6 2 2 1 2 1 0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 1 0];

        % budget of function evaluations and generations depending on dimension D
        MaxFES=D*20000;
        MaxGEN=D*2000;


        % input -- initial (fixed) strategy parameter setting
        input.budget            = MaxFES;
        input.maxIter           = MaxGEN;
        input.delta             = 10^-4;                    % error margin for equality constraints
        input.runs              = 25;                       % number of repetitions
        input.dim               = D;


        % CMA specific paprameters
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


        problem.constr_fun_name = 'CEC2018';                    % objective function class

        strategy='epsMAgES';                                    % algorithm variant
%         foldername = ['CEC2018_RS_' strategy];
%         efn = exist(foldername);
%         if efn ~= 7
%             mkdir(foldername);
%         end

        disp(['ES variant ' strategy ' --- 25 independent runs!'])

        for k=1:28                                              % on each CEC 2018 test function do
            func_num=k;                                         % test function number
            initial_flag = 0;

            eval(['input.lb=Xmin' int2str(func_num) ';']);
            eval(['input.ub=Xmax' int2str(func_num) ';' ]);

            problem.upper_bounds    = input.ub';
            problem.lower_bounds    = input.lb';


            for j=1:input.runs                                  % perform multiple runs on each test function
                eval(['[tab]=' strategy '(problem,input,func_num);']); % run epsMAgES
                FitT(j,:)=[tab(1,:) tab(2,:) tab(3,:)];
            end

            Tab=build_stats(FitT,input);                        % build statistics according to specification of the CEC2018 competition 
            Stats(k,:)=[func_num Tab(1,:) Tab(2,:) Tab(3,:)];
        end
              
        filename = strcat(['epsMAgES_RS_on_CEC2017_D' num2str(input.dim) '.txt']);
        fileID = fopen(filename,'w');             
        fprintf(fileID,'%14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\t %14.8g\n',Stats);
        fclose(fileID)
        
%         save([foldername '/' strategy '_RS_on_CEC2017_D' num2str(D) '.mat'],'input','Stats','-v7')
        
    end
