import pandas as pd

def confusion_raw(growth_file,essen_dic):
    out=[]
    fko_TF= open(growth_file,"r")
    fko_TF_lines=fko_TF.readlines()
    for line in fko_TF_lines[1:]:
        gen=line.split(",")[0]
        f=line.split(",")[1].rstrip("\n")
        essen_value=essen_dic[gen]
        out.append([gen, f, essen_value])

    return out

def confusion_metrics(V_P, V_N, F_P, F_N):
    total= V_P + V_N + F_P + F_N
    accuracy=(V_P+V_N)/total
    error_rate=1-accuracy
    sensitivity=V_P/(V_P+F_N)
    False_positive_rate=F_P/(V_N+F_P)
    True_positive_rate=1-False_positive_rate
    precision=V_P/(F_P+V_P)
    prevalence=(F_N+V_P)/total

    rounder = lambda num: round(num, 2)
 
    # print("accuracy", rounder(accuracy))
    # print("error_rate", rounder(error_rate))
    # print("sensitivity", rounder(sensitivity))
    # print("False_positive_rate", rounder(False_positive_rate))
    # print("True_positive_rate", rounder(True_positive_rate))
    # print("precision", rounder(precision))
    # print("prevalence", rounder(prevalence))

    return [rounder(accuracy), rounder(error_rate), rounder(sensitivity), 
    rounder(False_positive_rate), rounder(True_positive_rate), rounder(precision), rounder(prevalence)]

def confusion(prom_essen_file, essen_value, esse_threshold, growth_threshold, type_essen="griffin"):

    confusion_raw_data=confusion_raw(prom_essen_file,essen_value)
    V_P=[]
    F_N=[]
    F_P=[]
    V_N=[]

    for triada in confusion_raw_data:
        gen=triada[0]
        growth=float(triada[1])
        essen_value=triada[2]

        if type_essen == 'griffin':
            # Verdadero Positivo - predicts that it grows (not essential) and is correct.
            if essen_value>esse_threshold and growth>growth_threshold:
                V_P.append(gen)
            # Falso Negativo
            if essen_value>esse_threshold and growth<growth_threshold:
                F_N.append(gen)
            # Falso Positivo
            if essen_value<esse_threshold and growth>growth_threshold:
                F_P.append(gen)
            # Verdadero Negativo
            if essen_value<esse_threshold and growth<growth_threshold:
                V_N.append(gen)

        if type_essen == 'loerger':
            if (essen_value == "NE" or essen_value == "GA") and growth>growth_threshold:
                V_P.append(gen)
            if (essen_value == "NE" or essen_value == "GA") and growth<growth_threshold:
                F_N.append(gen)
            if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth>growth_threshold:
                F_P.append(gen)
            if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth<growth_threshold:
                V_N.append(gen)

    # print(type_essen)
    # print("VP",len(V_P))
    # print("VN",len(V_N))
    # print("FP",len(F_P))
    # print("FN",len(F_N))
    
    return confusion_metrics(len(V_P), len(V_N), len(F_P), len(F_N))


loerger_file = "mbo002173137st3.xlsx"
loerger_excel = pd.read_excel(loerger_file,skiprows = 1,keep_default_na=False)
loerger_finalCall = {loerger_excel.loc[idx, 'ORF ID']:  loerger_excel.loc[idx, 'Final Call'] for idx in range(loerger_excel.shape[0])}
loerger_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/f_DeJesus.txt"

griffin_file = "ppat.1002251.s002.xlsx"
griffin_excel = pd.read_excel(griffin_file,skiprows = 9,keep_default_na=False)
griffin_pvalue = {griffin_excel.loc[idx, 'Locus']:  griffin_excel.loc[idx, 'p value'] for idx in range(griffin_excel.shape[0])}
griffin_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/f_Griffin.txt"


# confusion(griffin_fko_TF_file, griffin_pvalue, esse_threshold=0.1, growth_threshold=0.5*0.0584)
# confusion( loerger_fko_TF_file, loerger_finalCall, esse_threshold=0.1, growth_threshold=0.5*0.0485, type_essen='loerger')

grafical_list=[]
#[round(x * 0.01, 2) for x in range(1, 100,1)]
#[round(x * 0.01, 1) for x in range(1, 100)]
intervals=[round(x * 0.01, 2) for x in range(1, 100,1)]
for x in intervals:
    #grafical_list.append(confusion(griffin_fko_TF_file, griffin_pvalue, esse_threshold=0.1, growth_threshold=x*0.0584))
    grafical_list.append(confusion( loerger_fko_TF_file, loerger_finalCall, esse_threshold=0.1, growth_threshold=x*0.0485, type_essen='loerger'))

accuracy_data=[]
error_rate_data=[]
sensitivity_data=[]
False_positive_rate_data=[]
True_positive_rate_data=[]
specificy_data=[] 
precision_data=[]
prevalence_data=[]
for element in grafical_list:
    accuracy_data.append(element[0])
    error_rate_data.append(element[1])
    sensitivity_data.append(element[2])
    False_positive_rate_data.append(element[3])
    True_positive_rate_data.append(element[4])
    precision_data.append(element[5])
    prevalence_data.append(element[6])

import matplotlib.pyplot as plt
plt.plot(
intervals, accuracy_data, 'r-', 
#intervals, error_rate_data, 'b-',
intervals, sensitivity_data, 'g-', 
intervals, False_positive_rate_data, 'b-', 
#intervals, True_positive_rate_data, 'b--',
intervals, precision_data, 'y-',
intervals, prevalence_data, 'k-',
)

plt.legend([
'accuracy_data', 
#'error_rate_data', 
'sensitivity_data', 
'False_positive_rate_data', 
#'True_positive_rate_data', 
'precision_data' , 
'prevalence_data'
])

# Datos paper iEK1011
#Grffin
#paper_iEK1011=confusion_metrics(579, 235, 161, 31)

#Loerger
paper_iEK1011=confusion_metrics(666, 221, 73, 45)


plt.plot(0.2, paper_iEK1011[0], "ro")
plt.plot(0.2, paper_iEK1011[2], "go")
plt.plot(0.2, paper_iEK1011[3], "bo")
plt.plot(0.2, paper_iEK1011[5], "yo")
plt.plot(0.2, paper_iEK1011[6], "ko")

plt.show()