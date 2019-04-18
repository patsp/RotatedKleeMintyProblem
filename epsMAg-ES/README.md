## epsMAg-ES
 Matlab code of the epsMAg-ES for constrained optimization problems

1. ### REFERENCE:  
   
  Hellwig, Michael; Beyer, Hans-Georg, "A Matrix Adaptation Evolution Strategy for Constrained Real-Parameter Optimization", Proceedings of IEEE Conference on Evolutionary Computation (CEC 2018), IEEE Xplore, 2018 (accepted).

2. ### NOTICE:  

  The epsMAg-ES code is made available for reproduction of reported results and testing. The use of the code is permitted subject to the condition of properly acknowledging this source (https://github.com/hellwigm/epsMAg-ES/) as well as citing the relevant papers.

3. ### The epsMAg-ES:  

  The epsMAg-ES represents a novel Evolution Strategy for constrained optimization that combines the the recently suggested Matrix Adaptation Evolution Strategy (MA-ES) with successful constraint handling techniques from the field of Differential Evolution. Being applied to the benchmark problems specified for the CEC 2018 competition on constrained single objective real-parameter optimization, the algorithm is able to find feasible solutions on more than 80% of the benchmark problems with high accuracy. 

4. ### Notice:
  This repository includes the benchmarking functions corresponding to the <b>"CEC 2018 competition on constrained single objective real-parameter optimization"</b> that are freely available on the website <href>http://www.ntu.edu.sg/home/epnsugan/index_files/CEC2018/CEC2018.htm</href>

5. ### Content:
* __Main_epsMAgES.m__ - Executable for running the epsMAg-ES on the correspondinng benchmarking problems
* __epsMAgES.m__ - Main component of the epsMAg-Es algorithm (based on the MA-ES)
* __eps_sort.m__ - Sorting routine for ranking candidate solutions w.r.t. the epsilon-level ordering
* __eps_rank.m__ - Subroutine of __eps_sort__
* __keep_range.m__ - Box-constraint handling method (reflection into the box)
* __gradientMutation.m__ - Implementation of the gradient-based repair step
* __build_stats.m__ - Postprocessing routine that builds the statistcs to be reported to the CEC 2018 competition
* __CEC2018.m__ - Constrained functions specified for the CEC 2017 and 2018 competitions
