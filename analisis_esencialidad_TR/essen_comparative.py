import pandas as pd

loerger_file = "mbo002173137st3.xlsx"
loerger_excel = pd.read_excel(loerger_file,skiprows = 1,keep_default_na=False)
loerger_finalCall = {loerger_excel.loc[idx, 'ORF ID']:  loerger_excel.loc[idx, 'Final Call'] for idx in range(loerger_excel.shape[0])}
loerger_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/f_DeJesus.txt"

griffin_file = "ppat.1002251.s002.xlsx"
griffin_excel = pd.read_excel(griffin_file,skiprows = 9,keep_default_na=False)
griffin_pvalue = {griffin_excel.loc[idx, 'Locus']:  griffin_excel.loc[idx, 'p value'] for idx in range(griffin_excel.shape[0])}
griffin_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/f_Griffin.txt"


def confusion_raw(growth_file,essen_dic):
    out=[]
    fko_TF= open(growth_file,"r")
    fko_TF_lines=fko_TF.readlines()
    for line in fko_TF_lines[1:]:
        gen=line.split(",")[0]
        f=line.split(",")[1].rstrip("\n")
        essen_value=essen_dic[gen]
        #print(gen,f,essen_value)
        out.append([gen, f, essen_value])
    return out

#Para griffin
griffin_confusion_raw=confusion_raw(griffin_fko_TF_file,griffin_pvalue)
# Segun paper iEK1011
griffin_threshold=0.1
# Segun paper iEK1011
growth_threshold=0.2*0.0584
V_P=[]
V_N=[]
F_P=[]
F_N=[]
for triada in griffin_confusion_raw:
    gen=triada[0]
    growth=float(triada[1])
    essen_value=triada[2]
    if essen_value>griffin_threshold and growth>growth_threshold:
        V_P.append(gen)
    if essen_value>griffin_threshold and growth<growth_threshold:
        V_N.append(gen)
    if essen_value<griffin_threshold and growth>growth_threshold:
        F_P.append(gen)
    if essen_value<griffin_threshold and growth<growth_threshold:
        F_N.append(gen)

print("Grffin")
print("VP",len(V_P))
print("NV",len(V_N))
print("FP",len(F_P))
print("FN",len(F_N))

#Para DeJesus/loerger
#Segun paper iEK1011
# Final call = NE o GA, Verdadero
# Final call = ES , ESD o GD, Negativo
loerger_confusion_raw=confusion_raw(loerger_fko_TF_file,loerger_finalCall)

growth_threshold=0.2*0.0485
V_P=[]
V_N=[]
F_P=[]
F_N=[]
for triada in loerger_confusion_raw:
    gen=triada[0]
    growth=float(triada[1])
    essen_value=triada[2]
    if (essen_value == "NE" or essen_value == "GA") and growth>growth_threshold:
        V_P.append(gen)
    if (essen_value == "NE" or essen_value == "GA") and growth<growth_threshold:
        V_N.append(gen)
    if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth>growth_threshold:
        F_P.append(gen)
    if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth<growth_threshold:
        F_N.append(gen)

print("Loerger")
print("VP",len(V_P))
print("NV",len(V_N))
print("FP",len(F_P))
print("FN",len(F_N))