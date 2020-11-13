
%Corro CobraToolBox
%initCobraToolbox(false);

% FBA - growth rate
%FBAsolution = optimizeCbModel(model,'max')

% Cargo los datos del paper de PROM
load('/home/agustin/FBA_Tesis/PROM_Chandrasekaran/mtbpromdata.mat');

fid = fopen('targets_filter.txt','r');
data = textscan(fid,'%s', 'Delimiter', '\n');
fclose(fid);
z_targets = data{1};

fid = fopen('regulator_filter.txt','r');
data = textscan(fid,'%s', 'Delimiter', '\n');
fclose(fid);
z_regulator = data{1};

z_litevidence = importfile_numeric('litevidence_filter.txt','r');
z_prob_prior = importfile_numeric('prob_prior_filter.txt','r');

clear ans fid data;


% Cargo el Modelo de BIGG, Soluciono problemade rules y rev.
load('iEK1008.mat')
rev=iEK1008.rev; 
clear ans iEK1008; 
iEK1008=readCbModel('iEK1008.mat'); 
iEK1008.rev=rev; 
clear rev;

% Chequeo Biomasa actual
checkObjective(iEK1008)
%printRxnFormula(iEK1008,'BIOMASS__2');

%% Cargo otro modelo, veo su reaccion de biomasa y la guardo
    load('/home/agustin/FBA_Tesis/PROM_trabajo/convertion/iEK1011_m7H10_media.mat');
    biomass = checkObjective(iEK1011);
    biomass_rxn_c = char(printRxnFormula(iEK1011,biomass));

% Cambio la biomasa
    %iEK1008 = addReaction(iEK1008, 'biomass_rxn_c','reactionFormula', biomass_rxn_c);
    %iEK1008 = changeObjective(iEK1008,'biomass_rxn_c');

% Chequeo el cambio de Biomasa
    checkObjective(iEK1008)
%%

addpath('/home/agustin/cobratoolbox/PROM_Chandrasekaran');

%PROM con datos del paper
%[f,v,status1,lostxns] = PROM(model, expression, expressionid, regulator, targets, litevidence,prob_prior);

%PROM con nuevos datos
%[f,v,status1,lostxns] = PROM(iEK1008, expression, expressionid, z_regulator, z_targets, z_litevidence, z_prob_prior);

% Corro PROM version 2
%function [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,v11,v12,KAPPA,DATATHRESHVAL,probtfgene,sizeflag)
[v11, v12] = fastFVA(iEK1008);
[f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(iEK1008,expression,expressionid,z_regulator,z_targets,z_litevidence,z_prob_prior,[],v11,v12,[],[],[],1)


%% Exporto variables
%Write the outputs in files
%fid = fopen('f.txt','w');
%CT = f;
%fprintf(fid,'%s\n', CT{:});
%fclose(fid);

%dlmwrite('f.txt', f, 'delimiter','\t','newline','pc','precision',13);

%% Analisis de resultados
diff_f= (f/0.0582)*100;
u_z_regulator=unique(z_regulator);
T = array2table(diff_f,'VariableNames',u_z_regulator);

YourArray = table2array(T);
Tt = array2table(YourArray.');
Tt.Properties.RowNames = T.Properties.VariableNames;
writetable(Tt,"reguladores_diff_f.txt",'WriteRowNames',true);