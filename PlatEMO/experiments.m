%   '-N'            <positive integer>  population size
%   '-M'            <positive integer>  number of objectives
%   '-D'            <positive integer>  number of variables
%	'-evaluation'   <positive integer>  maximum number of evaluations
%	'-algorithm'    <function handle>   algorithm function
%	'-problem'      <function handle>   problem function
%	'-operator'     <function handle>   operator function
%   '-mode'         <positive integer>  run mode (1.show result 2.save result 3.run outputFcn)
%   '-run'          <positive integer>  run No.
%   '-outputFcn'	<function handle>   function invoked after each generation when mode = 3
%	'-X_parameter'  <cell>              the parameter values of function X


% MyPar = parpool(10); % use 10 core to run the test
parfor i=1:30
    
    algorithm=@DNEA;
    runtest(i,algorithm);
    
    algorithm=@DNEAL;
    runtest(i,algorithm);
    
    algorithm=@MMOEAC;
    runtest(i,algorithm);
    
    algorithm=@HREA;
    runtest(i,algorithm);
    
    algorithm=@MMEAWI;
    runtest(i,algorithm);
    
    algorithm=@TriMOEATAR;
    runtest(i,algorithm);
    
    algorithm=@CPDEA;
    runtest(i,algorithm);
    
    algorithm=@MPMMEA;
    runtest(i,algorithm);
    
end
% delete(MyPar)

function runtest(runID,algorithm)
%% IDMP
for i=2:4
    for j=1:4
        popsize=100*i;
        evaluate=5000*i;
        command=['problem=@IDMPM' num2str(i) 'T' num2str(j) ';'];
        eval(command);
        main('-algorithm',algorithm,'-HREA_parameter',{0.2,0.5},'-problem',problem,'-N',popsize,'-evaluation',evaluate,'-run',runID,'-mode',2)
    end
end

%% IDMPe
for i=2:3
    for j=1:4
        popsize=100*i;
        evaluate=5000*i;
        command=['problem=@IDMPM' num2str(i) 'T' num2str(j) '_e;'];
        eval(command);
        main('-algorithm',algorithm,'-HREA_parameter',{0.5,0.5},'-problem',problem,'-N',popsize,'-evaluation',evaluate,'-run',runID,'-mode',2)
    end
end

%% MMF test
for j=1:8
    command=['problem=@MMF' num2str(j) ';'];
    eval(command);
    main('-algorithm',algorithm,'-problem',problem,'-N',200,'-evaluation',10000,'-run',runID,'-mode',2)
end

ep = [0.2,0.23,0.23,0.3,0.12,0.2,0.2];
for j=9:15
    if j>=13
        pop=300;eva=15000;
    else
        pop=200;eva=10000;
    end
    command=['problem=@MMF' num2str(j) ';'];
    eval(command);
    main('-algorithm',algorithm,'-HREA_parameter',{ep(j-8),0.5},'-problem',problem,'-N',pop,'-evaluation',eva,'-run',runID,'-mode',2)
end

%% Multi-Polygon
problem=@PolygonProblem;
count=0;
for i=[3,4,8,10]
    for j=[2,4,10,20,50]
        count=count+1;
        if i==3
            popsize=200;
        elseif i==4
            popsize=300;
        else
            popsize=400;
        end
        evaluate=max(50000,j*5000);
        main('-algorithm',algorithm,'-HREA_parameter',{0.2,0.5},'-problem',problem,'-M',i,'-N',popsize,'-D',j, '-evaluation',evaluate,'-run',j*1000+runID,'-mode',2)
     end
end

end