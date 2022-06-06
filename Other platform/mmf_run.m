% MyPar = parpool(10); % use 10 core to run the test
parfor i=1:30
    idrun(i);
end
% delete(MyPar)
% idrun(1);
function idrun(runID)
% runID=1;
%% Add path
addpath(genpath('MMF/'));
addpath(genpath('TruePSPFdata/'));
addpath(genpath('Indicator_calculation/'));
global fname
count=1;
for i=1:15
    if i >= 13%i,[13,14,15])
        popsize=300;%[200,210,210,156,210,230,135];
        evaluate=15000;%[10000,12000,15000,20000,20000,20000,20000];
    else
        popsize=200;%[200,210,210,156,210,230,135];
        evaluate=10000;%[10000,12000,15000,20000,20000,20000,20000];
    end
    
    count=count+1;
    fname=['MMF' num2str(i)];
    Max_evaluation=evaluate;
    Max_Gen=fix(Max_evaluation/popsize);
    
    var=feval(fname,'init',[],402);
    n_obj=var.M; M=2;
    n_var=var.D; D=2;
    xl=var.lower;
    xu=var.upper;
    
    fprintf('RUNID: %d Test function: %s Ä¿±ê£º%d Î¬¶È£º%d\n',runID, fname,n_obj,n_var);
    
    ffff=matfile(['Reference_PSPF_data\MMF' num2str(i) '_Reference_PSPF_data.mat']);
    PF=ffff.PF;
    PS=ffff.PS;
    
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
