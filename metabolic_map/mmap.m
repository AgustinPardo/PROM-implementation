initCobraToolbox(false) % false, as we don't want to update

% change to the directory of the tutorial
cd(fileparts(which('tutorial_ecoliCoreModel_part1.mlx')));

model = readCbModel('/home/agustin/FBA_Tesis/trabajo/iEK1008.mat');
mtb_coli_core = model; % Save the original model for later use

% Exporta un ".xls" con la data del modelo
%outmodel = writeCbModel(model, 'xls', 'core_model.xls');

map=readCbMap('/home/agustin/cobratoolbox/tutorials/reconstruction/ecoliCoreModel/ecoli_Textbook_ExportMap.txt');
%map=readCbMap('/home/agustin/FBA_Tesis/Paper_IEK1011/supplementaryMaterial/12918_2018_557_MOESM4_ESM/central_carbon.json')
options.zeroFluxWidth = 1;
options.rxnDirMultiplier = 10;
% Exporta un ".svg" denominado target
drawCbMap(map);

model = mtb_coli_core  ; % Starting with the original model
%model = changeRxnBounds(model,'EX_glc(e)',-10,'l'); % Set maximum glucose uptake
%model = changeRxnBounds(model,'EX_o2(e)',-30,'l'); % Set maximum oxygen uptake
%model = changeObjective(model,'Biomass_Ecoli_core_w_GAM'); % Set the objective function
FBAsolution = optimizeCbModel(model,'max') % FBA analysis

map=readCbMap('/home/agustin/cobratoolbox/tutorials/reconstruction/ecoliCoreModel/ecoli_Textbook_ExportMap.txt')
options.zeroFluxWidth = 0.1;
options.rxnDirMultiplier = 10;
drawFlux(map, model, FBAsolution.x, options); % Draw the flux values on the map "target.svg"

%rxns_flux=table(model.rxns,FBAsolution.v)
rxns_flux_json=jsonencode(containers.Map(model.rxns, FBAsolution.v))

fid = fopen('rxns_flux.json', 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, rxns_flux_json, 'char');
fclose(fid);
