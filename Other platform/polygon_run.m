% MyPar = parpool(10); % use 10 core to run the test
timm = tic;
parfor i=1:30
    idrun(1,tic);
end
% delete(MyPar)
% idrun(1);
function idrun(runID,startT)

%% Add path
addpath(genpath('PolygonTestProblem/'));
% addpath(genpath('TruePSPFdata/'));
addpath(genpath('Indicator_calculation/'));
global fname
count=1;
for i=[3,4,8,10]
    for j=[2,4,10,20,50]
        
        count=count+1;
        fname='PolygonProblem';
        if i==3
            popsize=200;
        elseif i==4
            popsize=300;
        else
            popsize=400;
        end
        Max_evaluation=max(50000,j*5000);
        Max_Gen=fix(Max_evaluation/popsize);
        
        Global.M=i;Global.D=j;
        var=feval(fname,'init',Global,402);
        n_obj=var.M; M=n_obj;
        n_var=var.D; D=n_var;
        xl=var.lower;
        xu=var.upper;
        
        fprintf('RUNID: %d Test function: %s 目标�?%d 维度�?%d 用时:%.2f\n',runID, fname,i,j,toc(startT));
        
        Global.N=10000;
        [PS,PF]=feval(fname,'PF',Global,402);
        PS=repmat(PS, 1, j/2);
        
        %% Search the PSs using DN
        t1=tic;
        disp('dnnsga2')
        [ps,pf]=DN_NSGAII(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        runtime=toc(t1);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\DN_NSGAII_' fname '_' num2str(M) '_' num2str(D) '-' num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf','runtime');
        
        %% Search the PSs using OMNI-opt
        t1=tic;
        disp('omniopt')
        [ps,pf]=Omni_Opt(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        runtime=toc(t1);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\Omni_Opt_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf','runtime');
        
        %% Search the PSs using MO_Ring_PSO_SCD
        t1=tic;
        disp('moringpsoscd')
        [ps,pf]=MO_Ring_PSO_SCD(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        runtime=toc(t1);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\MO_Ring_PSO_SCD_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf','runtime');
        
        %% Search the PSs using MO_PSO_MM
        t1=tic;
        disp('mopsomm')
        [ps,pf]=MO_PSO_MM(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        runtime=toc(t1);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\MO_PSO_MM_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf','runtime');
        
    end
end
end
