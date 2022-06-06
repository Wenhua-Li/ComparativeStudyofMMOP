% MyPar = parpool(10); % use 10 core to run the test
parfor i=1:30
    idrun(i);
end
% delete(MyPar)
% idrun(1);
function idrun(runID)
% runID=1;
%% Add path
addpath(genpath('IDMP_e/'));
addpath(genpath('TruePSPFdata/'));
addpath(genpath('Indicator_calculation/'));
global fname
count=1;
for i=2:3
    for j=1:4
        count=count+1;
        fname=['IDMPM' num2str(i) 'T' num2str(j) '_e'];
        popsize=i*100;
        Max_evaluation=i*5000;
        Max_Gen=fix(Max_evaluation/popsize);
        
        var.M=i;
        M=i;
        D=i;
        var=feval(fname,'init',var,402);
        var=[];
        var.M=i;
        n_obj=i;
        n_var=i;
        xl=-1.*ones(1,i);
        xu=1.*ones(1,i);
        
        fprintf('RUNID: %d Test function: %s Ä¿±ê£º%d Î¬¶È£º%d\n',runID, fname,i,i);
        
        ffff=matfile(['TruePSPFdata\IDMPM' num2str(i) 'T' num2str(j) '_e.mat']);
        PF=ffff.PFs;
        PS=ffff.PSs;
        
        %% Search the PSs using DN
        disp('dnnsga2')
        [ps,pf]=DN_NSGAII(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\DN_NSGAII_' fname '_' num2str(M) '_' num2str(D) '-' num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf');
        
         %% Search the PSs using OMNI-opt
         disp('omniopt')
        [ps,pf]=Omni_Opt(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\Omni_Opt_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf');
        
       %% Search the PSs using MO_Ring_PSO_SCD
       disp('moringpsoscd')
        [ps,pf]=MO_Ring_PSO_SCD(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\MO_Ring_PSO_SCD_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf');
        
        %% Search the PSs using MO_PSO_MM
        disp('mopsomm')
        Max_Gen=fix(Max_evaluation/popsize);
        [ps,pf]=MO_PSO_MM(fname,xl,xu,n_obj,popsize,Max_Gen,var);
        IGD=IGD_calculation(pf,PF);
        IGDx=IGD_calculation(ps,PS);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;% Eq. (8) in the paper
        r_file=['data\MO_PSO_MM_' fname '_' num2str(M) '_' num2str(D) '-'  num2str(runID) '.mat'];
        save(r_file,'IGD','IGDx','CR','PSP','ps','pf');

    end
end
end
