% MyPar = parpool(10); % use 10 core to run the test
parfor i=1:30
    idrun(i);
end
% delete(MyPar)
% idrun(1);
function idrun(runID)
% runID=1;
%% Add path
addpath(genpath('IDMP/'));
addpath(genpath('TruePSPFdata/'));
addpath(genpath('Indicator_calculation/'));
% popsizee=[200,210,210,156,210,230,135];
% evaluate=[10000,12000,15000,20000,20000,20000,20000];
global fname
count=1;
for i=2:4
    for j=1:4
        count=count+1;
        fname=['IDMPM' num2str(i) 'T' num2str(j)];
        popsize=i*100;
        Max_evaluation=i*5000;
        Max_Gen=fix(Max_evaluation/popsize);
        
        var=feval(fname,'init',[],402);
        n_obj=var.M; M=i;
        n_var=var.D; D=i;
        xl=var.lower;
        xu=var.upper;
        
        fprintf('RUNID: %d Test function: %s Ä¿±ê£º%d Î¬¶È£º%d\n',runID, fname,i,i);
        
        ffff=matfile(['TruePSPFdata\IDMPM' num2str(i) 'T' num2str(j) '.mat']);
        PF=ffff.PF;
        PS=ffff.PSS;
        
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
