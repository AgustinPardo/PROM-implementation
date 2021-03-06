
%Corro CobraToolBox
%initCobraToolbox(false);

% FBA - growth rate
%FBAsolution = optimizeCbModel(model,'max')

% Cargo los datos del paper de PROM
load('/home/agustin/FBA_Tesis/PROM_Chandrasekaran/mtbpromdata.mat');

% Cargo los targets
fid = fopen('/home/agustin/FBA_Tesis/PROM_trabajo/interseccion_reg_met/ernesto_vs_iEK1011/target_filter.txt','r');
data = textscan(fid,'%s', 'Delimiter', '\n');
fclose(fid);
z_targets = data{1};
clearvars fid data;

% Cargo los reguladores
fid = fopen('/home/agustin/FBA_Tesis/PROM_trabajo/interseccion_reg_met/ernesto_vs_iEK1011/regulator_filter.txt','r');
data = textscan(fid,'%s', 'Delimiter', '\n');
fclose(fid);
z_regulator = data{1};
clearvars fid data;

% Cargo z_litevidence
filename = '/home/agustin/FBA_Tesis/PROM_trabajo/interseccion_reg_met/ernesto_vs_iEK1011/litevidence_filter.txt';
delimiter = {''};
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
z_litevidence = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

% Cargo z_prob_prior
filename = '/home/agustin/FBA_Tesis/PROM_trabajo/interseccion_reg_met/ernesto_vs_iEK1011/prob_prior_filter.txt';
delimiter = {''};
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
z_prob_prior = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans data;

% Nueva matriz de expresion
expression_colombos_1021=dlmread('/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/expression_colombos_1021.txt')
% Nuevo expressionid
fid = fopen('/home/agustin/FBA_Tesis/PROM_trabajo/COLOMBOS_data/expressionid_colombos_1021.txt','r');
data = textscan(fid,'%s', 'Delimiter', '\n');
fclose(fid);
expressionid_colombos_1021 = data{1};
clear ans fid data;

% Cargo el Modelo. Soluciono problemade rules y rev.
% Al cargarlo con load, se soluciona rules y rev se cargan. Si lo cargo con
% readCbodel no se soluciona
load('/home/agustin/FBA_Tesis/PROM_trabajo/convertion/iEK1011_inVivo_media.mat')
rev=iEK1011.rev;
clear ans iEK1011; 
iEK1011=readCbModel('/home/agustin/FBA_Tesis/PROM_trabajo/convertion/iEK1011_inVivo_media.mat'); 
iEK1011.rev=rev; 
clear rev;
% Chequeo Biomasa actual
%checkObjective(iEK1011)
%printRxnFormula(iEK1008,'BIOMASS__2');
FBAsolution = optimizeCbModel(iEK1011,'max')


% Cargo otro modelo, veo su reaccion de biomasa y la guardo
    load('/home/agustin/FBA_Tesis/PROM_trabajo/convertion/iEK1011_inVivo_media.mat');
    biomass = checkObjective(iEK1011);
    biomass_rxn_c = char(printRxnFormula(iEK1011,biomass));

Cambio la biomasa
    iEK1011 = addReaction(iEK1011, 'biomass_rxn_c','reactionFormula', biomass_rxn_c);
    iEK1011 = changeObjective(iEK1011,'biomass_rxn_c');

Chequeo el cambio de Biomasa
    checkObjective(iEK1008)
%%

% Veo el medio. Para eso tengo que ver las reacciones de exchange EX
printConstraints(iEK1011,-500, +500)
%%

addpath('/home/agustin/cobratoolbox/PROM_Chandrasekaran');

%PROM con datos del paper
%[f,v,status1,lostxns] = PROM(model, expression, expressionid, regulator, targets, litevidence,prob_prior);

%PROM con nuevos datos
%[f,v,status1,lostxns] = PROM(iEK1008, expression, expressionid, z_regulator, z_targets, z_litevidence, z_prob_prior);

%% Corro PROM version 2
%function [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,v11,v12,KAPPA,DATATHRESHVAL,probtfgene,sizeflag)

[v11, v12] = fastFVA(iEK1011);
[f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(iEK1011,expression,expressionid,z_regulator, z_targets,z_litevidence,z_prob_prior,[],v11,v12,[],[],[],1);

%% PROM con un for, guardando las salidas
u_z_regulator=unique(z_regulator);
f_table=table(u_z_regulator);
percentil_tresh=[0,0.01,0.05,0.1,0.2,0.25,0.33,0.4,0.5,0.75,1]
for i = percentil_tresh
    [v11, v12] = fastFVA(iEK1011);
    [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(iEK1011,expression,expressionid,z_regulator,z_targets,z_litevidence,z_prob_prior,[],v11,v12,[],[i],[],1);
    T=table(transpose(f),'VariableNames', "t"+string(i));
    f_table=[f_table T];
end
%% Exporto las probabilidades (probtfgene). Tengo que unir z_regulator, z_targets y probtfgene
T = table(z_regulator, z_targets, probtfgene);
writetable(T,"probtfgene_Ernesto_colombos.txt");

%% Exporto variables
%Write the outputs in files
fid = fopen('iEK1011_genes.txt','w');
CT = iEK1011.genes;
fprintf(fid,'%s\n', CT{:});
fclose(fid);

dlmwrite('f.txt', f, 'delimiter','\t','newline','pc','precision',13);

% Para guardar una matriz
% dlmwrite('file.txt',expression)

%Export matriz como tabla
% exp_table = array2table(expression(1:4295, :));
% exp_table.Properties.RowNames=exp_id;
% writetable(exp_table,"exp_PROM_Original.txt",'WriteRowNames',true);

%% Analisis de resultados
% diff
FBAsolution = optimizeCbModel(iEK1011,'max')
diff_f= (f/FBAsolution.f)*100;
regulator=unique(z_regulator);
T = array2table(diff_f,'VariableNames',regulator);

YourArray = table2array(T);
Tt = array2table(YourArray.');
Tt.Properties.RowNames = T.Properties.VariableNames;
writetable(Tt,"ei437_invivo.txt",'WriteRowNames',true);

% f solo
t=table(regulator,transpose(f));
writetable(t,"ei437_invivol.txt",'WriteRowNames',true);
