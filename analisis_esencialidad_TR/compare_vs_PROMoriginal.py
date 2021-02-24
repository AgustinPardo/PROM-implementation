import pandas
candidate_essential_PROM=["Rv1395","Rv3291c","Rv2069","Rv3223c","Rv1657","Rv3575c"]
evidence_esential_PROM_correct=["Rv0001","Rv0485","Rv3676","Rv3414c","Rv2711"]
evidence_esential_PROM_incorrect=["Rv1027c"]
non_esential_PROM_correct=["Rv0117","Rv0212c","Rv0353","Rv0491","Rv0586",
                            "Rv0844c","Rv1221","Rv1909c","Rv1931c","Rv2359","Rv2720",
                            "Rv3080c","Rv3133c","Rv3279c","Rv3286c","Rv3574","Rv1785c",
                            "Rv1267c"]
all_PROM_analized_genes=candidate_essential_PROM+evidence_esential_PROM_correct+evidence_esential_PROM_incorrect+non_esential_PROM_correct
PROM_essential_genes=candidate_essential_PROM+evidence_esential_PROM_correct

table = {'Type of prediction': ['candidate_essential_PROM','evidence_esential_PROM_correct',
                                'evidence_esential_PROM_incorrect','non_esential_PROM_correct', "Totales", 'Total escenciales', 'Nuevos escenciales'],
        }
df_table = pandas.DataFrame(table, columns = ['Type of prediction'])

files_list=["/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/diff_f_DeJesus_ei437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/diff_f_Griffin_ei437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_437/diff_f_m7H10_ei437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/diff_f_DeJesus_eic.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/diff_f_Griffin_eic.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Ernesto_iEK1011_colombos/diff_f_m7H10_eic.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_437/diff_f_DeJesus_si437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_437/diff_f_Griffin_si437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_437/diff_f_m7H10_si437.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_colombos/diff_f_DeJesus_sic.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_colombos/diff_f_Griffin_sic.txt",
            "/home/agustin/FBA_Tesis/PROM_trabajo/analisis_esencialidad_TR/Sanz_iEK1011_colombos/diff_f_m7H10_sic.txt"
            ]
for file in files_list:
    label= file.split("/")[-1].split(".")[0]
    df2 = pandas.read_csv(file)
    dictionary=pandas.Series(df2.Var1.values,index=df2.Row).to_dict()
    #95
    growth_tresh=99.5

    candidate_essential_PROM_P=[]
    candidate_essential_PROM_N=[]
    evidence_esential_PROM_correct_P=[]
    evidence_esential_PROM_correct_N=[]
    evidence_esential_PROM_incorrect_P=[]
    evidence_esential_PROM_incorrect_N=[]
    non_esential_PROM_correct_P=[]
    non_esential_PROM_correct_N=[]

    for element in candidate_essential_PROM:
        prediction_value=float(dictionary[element])
        if prediction_value < growth_tresh:
            candidate_essential_PROM_P.append(element)
        else:
            candidate_essential_PROM_N.append(element)

    for element in evidence_esential_PROM_correct:
        prediction_value=float(dictionary[element])
        if prediction_value < growth_tresh:
            evidence_esential_PROM_correct_P.append(element)
        else:
            evidence_esential_PROM_correct_N.append(element)

    for element in evidence_esential_PROM_incorrect:
        prediction_value=float(dictionary[element])
        if prediction_value < growth_tresh:
            evidence_esential_PROM_incorrect_P.append(element)
        else:
            evidence_esential_PROM_incorrect_N.append(element)

    for element in non_esential_PROM_correct:
        prediction_value=float(dictionary[element])
        if prediction_value > growth_tresh:
            non_esential_PROM_correct_P.append(element)
        else:
            non_esential_PROM_correct_N.append(element)

    # print(label)
    # print("candidate_essential_PROM")
    # print("candidate_essential_PROM P: " + str(len(candidate_essential_PROM_P)))
    # print("candidate_essential_PROM N: " + str(len(candidate_essential_PROM_N)))

    # print("evidence_esential_PROM_correct")
    # print("evidence_esential_PROM_correct P: " + str(len(evidence_esential_PROM_correct_P)))
    # print("evidence_esential_PROM_correct N: " + str(len(evidence_esential_PROM_correct_N)))

    # print("evidence_esential_PROM_incorrect")
    # print("evidence_esential_PROM_incorrect P: " + str(len(evidence_esential_PROM_incorrect_P)))
    # print("evidence_esential_PROM_incorrect N: " + str(len(evidence_esential_PROM_incorrect_N)))

    # print("non_esential_PROM_correct")
    # print("non_esential_PROM_correct P: " + str(len(non_esential_PROM_correct_P)))
    # print("non_esential_PROM_correct N: " + str(len(non_esential_PROM_correct_N)))

    positivos = len(candidate_essential_PROM_P)+ len(evidence_esential_PROM_correct_P) + len(evidence_esential_PROM_incorrect_P) + len(non_esential_PROM_correct_P)
    negativos = len(candidate_essential_PROM_N)+ len(evidence_esential_PROM_correct_N) + len(evidence_esential_PROM_incorrect_N) + len(non_esential_PROM_correct_N)
    
    # print("Resumen")
    # print("Positivos: " + str(positivos))
    # print("Negativos: " + str(negativos))

    total_esenciales=[]
    new_esenciales=[]
    for key,value in dictionary.items():
        if value < growth_tresh:
            total_esenciales.append(key)
            if key not in PROM_essential_genes:
                new_esenciales.append(key)
    # Teniendo en cuenta la predicion sin compararla con
    # print("Total escenciales: "+ str(len(total_esenciales)))
    # print("Nuevos esenciales: "+ str(len(new_esenciales)))
    print(positivos)
    data = [len(candidate_essential_PROM_P), len(evidence_esential_PROM_correct_P), 
            len(evidence_esential_PROM_incorrect_P), len(non_esential_PROM_correct_P), positivos, len(total_esenciales), len(new_esenciales) ]
    df_table[label+" Positivo"] = data

    data = [len(candidate_essential_PROM_N), len(evidence_esential_PROM_correct_N), 
            len(evidence_esential_PROM_incorrect_N), len(non_esential_PROM_correct_N), negativos,"-", "-"]
    df_table[label+" Negativo"] = data

print(df_table)
df_table= df_table.T
df_table.to_csv('compare_vs_PROMoriginal.csv')


