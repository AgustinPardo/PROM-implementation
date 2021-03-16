import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

def confusion_raw(growth_file,essen_dic):
    out=[]
    fko_TF= open(growth_file,"r")
    fko_TF_lines=fko_TF.readlines()
    rounder = lambda num: round(num, 4)
    for line in fko_TF_lines[1:]:
        gen=line.split(",")[0]
        f=rounder(float(line.split(",")[1].rstrip("\n")))
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
    try:
        precision=V_P/(F_P+V_P)
    except ZeroDivisionError:
        precision=0
    prevalence=(F_N+V_P)/total

    rounder = lambda num: round(num, 2)
 
    print("accuracy", rounder(accuracy))
    print("error_rate", rounder(error_rate))
    print("sensitivity", rounder(sensitivity))
    print("False_positive_rate", rounder(False_positive_rate))
    print("True_positive_rate", rounder(True_positive_rate))
    print("precision", rounder(precision))
    print("prevalence", rounder(prevalence))

    return [rounder(accuracy), rounder(error_rate), rounder(sensitivity), 
            rounder(False_positive_rate), rounder(True_positive_rate), rounder(precision), 
            rounder(prevalence),V_P/total, V_N/total, F_P/total, F_N/total]

FN_conjuto=[]

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

        # Filtro genes que me interesan
        if gen in core_genes_uniprot_mycobrowser:

            if type_essen == 'griffin':
                # Verdadero Negativo - predicts that it grows (not essential) and is correct.
                # print(growth, growth_threshold)
                if essen_value>esse_threshold and growth>growth_threshold:
                    V_N.append(gen)
                # Falso Positivo
                if essen_value>=esse_threshold and growth<=growth_threshold:
                    F_P.append(gen)
                # Falso Negativo
                if essen_value<esse_threshold and growth>growth_threshold:
                    F_N.append(gen)
                # Verdadero Positivo
                if essen_value<=esse_threshold and growth<=growth_threshold:
                    V_P.append(gen)

        if type_essen == 'loerger':
#Final Call key: Essential (ES), Essential Domain (ESD), Growth-Defect (GD), Non-Essential (NE),  Growth-Advantage (GA), and Uncertain (for short empty genes).
            if (essen_value == "NE" or essen_value == "GA") and growth>growth_threshold:
                V_P.append(gen)
            if (essen_value == "NE" or essen_value == "GA") and growth<growth_threshold:
                F_N.append(gen)
            if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth>growth_threshold:
                F_P.append(gen)
            if (essen_value == "ES" or essen_value == "ESD" or essen_value == "GD") and growth<growth_threshold:
                V_N.append(gen)

            
    print(prom_essen_file.split("/")[-1])
    print(type_essen)
    if type_essen == "griffin":
        print(growth_threshold)
    print("VP",len(V_P))
    print("VN",len(V_N))
    print("FP",len(F_P))
    print("FN",len(F_N))

    print("VP", V_P)
    print("VN", V_N)
    print("FP", F_P)
    print("FN", F_N)
    

    return confusion_metrics(len(V_P), len(V_N), len(F_P), len(F_N))


loerger_file = "mbo002173137st3.xlsx"
loerger_excel = pd.read_excel(loerger_file,skiprows = 1,keep_default_na=False)
loerger_finalCall = {loerger_excel.loc[idx, 'ORF ID']:  loerger_excel.loc[idx, 'Final Call'] for idx in range(loerger_excel.shape[0])}
loerger_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/f_Griffin_ei437.txt"

griffin_file = "ppat.1002251.s002.xlsx"
griffin_excel = pd.read_excel(griffin_file,skiprows = 9,keep_default_na=False)
griffin_pvalue = {griffin_excel.loc[idx, 'Locus']:  griffin_excel.loc[idx, 'p value'] for idx in range(griffin_excel.shape[0])}
griffin_fko_TF_file="/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/f_Griffin_eic.txt"

# griffin_sassetti_data= {griffin_excel.loc[idx, 'Locus']:  griffin_excel.loc[idx, '(Sassetti et al 2003)'] for idx in range(griffin_excel.shape[0])}
# fko_TF= open(griffin_fko_TF_file,"r")
# fko_TF_lines=fko_TF.readlines()[1:]
# for line in fko_TF_lines:
#     gen=line.split(",")[0]
#     data=griffin_sassetti_data[gen]
#     print(gen, data, sep = " ")

#confusion( griffin_fko_TF_file, griffin_pvalue, esse_threshold=0.1, growth_threshold=0.95*0.0584)


# confusion( loerger_fko_TF_file, loerger_finalCall, esse_threshold=0.1, growth_threshold=0.2*0.0485, type_essen='loerger')

loerger_files=["/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/f_DeJesus_ei437.txt",
                "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/f_DeJesus_eic.txt",
                "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_437/f_DeJesus_si437.txt",
                "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_colombos/f_DeJesus_sic.txt"]

griffin_files=["/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/f_Griffin_ei437.txt",
                # "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/f_Griffin_eic.txt",
                # "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_437/f_Griffin_si437.txt",
                # "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_colombos/f_Griffin_sic.txt"
                ]

core_genes=["Rv3246c","Rv0001","Rv1221","Rv1475c","Rv2711","Rv0212c",
            "Rv0494","Rv0981","Rv1395","Rv1657","Rv1846c","Rv2213",
            "Rv2710","Rv3286c","Rv3291c","Rv3575c","Rv0117","Rv0182c","Rv0353",
            "Rv0445c","Rv0465c","Rv0576","Rv0586","Rv0602","Rv0792c","Rv0818",
            "Rv0823c","Rv0844c","Rv1027c","Rv1033c","Rv1267c","Rv1287","Rv1316c","Rv1359",
            "Rv1675c","Rv1785c","Rv1956","Rv1963c","Rv1985c","Rv1994c","Rv2021c","Rv2034",
            "Rv2175c","Rv2704","Rv2718c","Rv2788","Rv3080c","Rv3133c","Rv3164c","Rv3328c","Rv3334",
            "Rv3414c","Rv3416","Rv3557c","Rv3678c","Rv3681c","Rv3692","Rv3744"]

core_genes_uniprot=["Rv0981","Rv2710","Rv3286c","Rv3223c","Rv0182c","Rv0844c","Rv1027c",
                    "Rv1956","Rv1963c","Rv3133c","Rv3414c","Rv3416","Rv0465c","Rv3080c","Rv3328c","Rv3681c",
                    "Rv2069","Rv3911","Rv3676","Rv3003c","Rv3246c","Rv0485","Rv1931c","Rv2745c","Rv3082c",
                    "Rv3849","Rv1221"]

core_genes_uniprot_mycobrowser= ["Rv1221","Rv2711","Rv2374c","Rv0491","Rv3246c","Rv0485","Rv1931c","Rv2745c",
                            "Rv3082c","Rv3849","Rv0182c","Rv0844c","Rv1027c","Rv0981","Rv1956","Rv1963c","Rv3133c","Rv3414c","Rv3416",
                            "Rv2069","Rv3911","Rv2359","Rv1909c","Rv1316c","Rv1675c","Rv1994c","Rv2718c","Rv0465c","Rv3080c","Rv3328c","Rv3681c",
                            "Rv2710","Rv3286c","Rv3223c","Rv0576","Rv0212c"]

# CAMBIAR el set de archivos a usar
for file in griffin_files:
    grafical_list=[]
    #[round(x * 0.01, 2) for x in range(1, 100,1)]
    #[round(x * 0.01, 1) for x in range(1, 100)]
    intervals=[round(x * 0.01, 2) for x in range(0, 101,1)]

    for x in intervals:
        label= file.split("/")[-1].split(".")[0]
        # CAMBIAR!
        grafical_list.append(confusion(file, griffin_pvalue, esse_threshold=0.1, growth_threshold=x*0.0584))
        #grafical_list.append(confusion( file, loerger_finalCall, esse_threshold=0.1, growth_threshold=x*0.0485, type_essen='loerger'))

    accuracy_data=[]
    error_rate_data=[]
    sensitivity_data=[]
    False_positive_rate_data=[]
    True_positive_rate_data=[]
    specificy_data=[] 
    precision_data=[]
    prevalence_data=[]
    V_P=[]
    V_N=[]
    F_P=[]
    F_N=[]
    for element in grafical_list:
        accuracy_data.append(element[0])
        error_rate_data.append(element[1])
        sensitivity_data.append(element[2])
        False_positive_rate_data.append(element[3])
        True_positive_rate_data.append(element[4])
        precision_data.append(element[5])
        prevalence_data.append(element[6])
        V_P.append(element[7])
        V_N.append(element[8])
        F_P.append(element[9])
        F_N.append(element[10])

    # ROC curve
    #print('AUC: {}'.format(auc(False_positive_rate_data, sensitivity_data)))
    # plt.figure(figsize=(10, 8))
    # lw = 2
    # plt.plot(False_positive_rate_data, sensitivity_data, color='darkorange',
    #          lw=lw, label='ROC curve')
    # plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.yticks([i/20.0 for i in range(21)])
    # plt.xticks([i/20.0 for i in range(21)])
    # plt.xlabel('False Positive Rate')
    # plt.ylabel('True Positive Rate')
    # plt.title(label)
    # plt.plot([], [], ' ', label='AUC: {}'.format(auc(False_positive_rate_data, sensitivity_data)))
    # plt.legend(loc='lower right')
    # plt.savefig("ROC_"+label+".png", dpi = 200)
    # plt.clf()
    

# #     #Curvas de metricas de confusion cambiando el growth rate threshold
# #     plt.plot(
# #     intervals, accuracy_data, 'r-', 
# #     #intervals, error_rate_data, 'b-',
# #     intervals, sensitivity_data, 'g-', 
# #     intervals, False_positive_rate_data, 'b-', 
# #     #intervals, True_positive_rate_data, 'b--',
# #     intervals, precision_data, 'y-',
# #     #intervals, prevalence_data, 'k-',
# #     )

# #     plt.title(label)

# #     plt.ylim([0, 1])

# #     plt.xlabel("Grow rate threshold")
# #     plt.ylabel("Prediction metric")

# #     plt.legend([
# #     'accuracy_data', 
# #     #'error_rate_data', 
# #     'sensitivity_data', 
# #     'False_positive_rate_data', 
# #     #'True_positive_rate_data', 
# #     'precision_data' , 
# #     #'prevalence_data'
# #     ])

# #     # #Datos paper iEK1011
# #     # #Griffin
# #     # paper_iEK1011=confusion_metrics(579, 235, 161, 31)

# #     # #Loerger
# #     # paper_iEK1011=confusion_metrics(666, 221, 73, 45)

# #     # plt.plot(0.2, paper_iEK1011[0], "ro")
# #     # plt.plot(0.2, paper_iEK1011[2], "go")
# #     # plt.plot(0.2, paper_iEK1011[3], "bo")
# #     # plt.plot(0.2, paper_iEK1011[5], "yo")
# #     # plt.plot(0.2, paper_iEK1011[6], "ko")

# #     #plt.savefig(label+"VPng_ng"+".png", dpi = 200)
# #     plt.clf()

# # #     # #ROC curve
# # #     # plt.plot( False_positive_rate_data, sensitivity_data )
# # #     # plt.plot([-1, 1], [-1, 1], ls="--", c=".3")
# # #     # plt.title(label)
# # #     # plt.ylim([0, 1])
# # #     # plt.xlim([0, 1])
# # #     # plt.xlabel("False positive rate")
# # #     # plt.ylabel("Sensitivity (True positive rate)")

# # #     # plt.savefig("ROC_"+label+".png", dpi = 200)
# # #     # plt.clf()

    # Curvas de V_P, V_N, F_P, F_N sobre total
    # Curvas de metricas de confusion cambiando el growth rate threshold
    plt.plot(
    intervals, V_P, 'r-', 
    intervals, V_N, 'g-', 
    intervals, F_P, 'b-', 
    intervals, F_N, 'y-',
    )

    plt.title(label)

    plt.ylim([0, 1])

    plt.xlabel("Grow rate threshold")
    plt.ylabel("Prediction")

    plt.legend([
    'V_P', 
    'V_N', 
    'F_P',
    'F_N' ,
    ])
    plt.xticks([i/20.0 for i in range(21)])
    plt.xticks(fontsize=8, rotation=90)
    plt.grid()
    plt.savefig(label+"_core_uniprot_mycobrowser"+".png", dpi = 200)
    plt.clf()



griffin_threshold_by_graphic=[0.85,0.85,0.65,0.8]
griffin_threshold_by_graphic=[0.65]
DeJesus_threshold_by_graphic=[0.8,0.75,0.525,0.675]
# Calculo la mtriz de confusion para el threshold seleccionado y veo quienes son los FP
for index in range(len(griffin_files)):
    confusion(loerger_files[index], griffin_pvalue, esse_threshold=0.1, growth_threshold=griffin_threshold_by_graphic[index]*0.0584)