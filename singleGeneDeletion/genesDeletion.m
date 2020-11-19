iEK1011 = readCbModel('/home/agustin/FBA_Tesis/PROM_trabajo/convertion/iEK1011_deJesusEssen_media.mat');
[grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion(iEK1011);

genes=iEK1011.genes;
t=table(genes,grRateKO);
writetable(t,"DeJesus_Essen_metabolic.txt",'WriteRowNames',true);
