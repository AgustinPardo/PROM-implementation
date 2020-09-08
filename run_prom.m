clear;

load('/home/agustin/cobratoolbox/PROM_Chandrasekaran/mtbpromdata.mat');

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

load('iEK1008.mat');

rev=iEK1008.rev;

clear ans iEK1008;

iEK1008=readCbModel('iEK1008.mat');

iEK1008.rev=rev;

clear rev;

addpath('/home/agustin/cobratoolbox/PROM_Chandrasekaran');

%PROM con datos del paper
%[f,v,status1,lostxns] = PROM(model, expression, expressionid, regulator, targets, litevidence,prob_prior);

%PROM con nuevos datos
[f,v,status1,lostxns] = PROM(iEK1008, expression, expressionid, z_regulator, z_targets, z_litevidence, z_prob_prior);


%fid = fopen('expressionid.txt','w');
%CT = expressionid;
%fprintf(fid,'%s\n', CT{:});
%fclose(fid);
