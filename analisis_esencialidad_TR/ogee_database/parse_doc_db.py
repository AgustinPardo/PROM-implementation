import pandas
ernesto_FT_genes=["Rv0001","Rv0117","Rv0182c","Rv0212c","Rv0260c","Rv0348","Rv0353",
                "Rv0445c","Rv0465c","Rv0485","Rv0491","Rv0494","Rv0576","Rv0586","Rv0602c","Rv0757",
                "Rv0792c","Rv0818","Rv0823c","Rv0844c","Rv0981","Rv1027c","Rv1033c","Rv1221","Rv1267c",
                "Rv1287","Rv1316c","Rv1359","Rv1379","Rv1395","Rv1475c","Rv1657","Rv1675c","Rv1785c","Rv1846c",
                "Rv1909c","Rv1931c","Rv1956","Rv1963c","Rv1985c","Rv1994c","Rv2017","Rv2021c","Rv2034","Rv2069",
                "Rv2166c","Rv2175c","Rv2213","Rv2359","Rv2374c","Rv2704","Rv2710","Rv2711","Rv2718c","Rv2720","Rv2745c",
                "Rv2788","Rv3002c","Rv3003c","Rv3080c","Rv3082c","Rv3133c","Rv3164c","Rv3223c","Rv3246c","Rv3279c","Rv3286c",
                "Rv3291c","Rv3328c","Rv3334","Rv3414c","Rv3416","Rv3557c","Rv3574","Rv3575c","Rv3676","Rv3678c","Rv3681c","Rv3692",
                "Rv3744","Rv3849","Rv3911"
                ]

df = pandas.read_csv("Mycobacterium tuberculosis H37Rv_genes.csv")
sassetti2003Invivo = df[df['dataset'] == 168]
sassetti2003Invivo_dic=pandas.Series(sassetti2003Invivo.essentiality.values,index=sassetti2003Invivo.locus).to_dict()
Lamichhane2003 = df[df['dataset'] == 169]
Lamichhane2003_dic=pandas.Series(Lamichhane2003.essentiality.values,index=Lamichhane2003.locus).to_dict()
sassetti2003Invitro = df[df['dataset'] == 170]
sassetti2003Invitro_dic=pandas.Series(sassetti2003Invitro.essentiality.values,index=sassetti2003Invitro.locus).to_dict()
Rubin2012 = df[df['dataset'] == 171]
Rubin2012_dic=pandas.Series(Rubin2012.essentiality.values,index=Rubin2012.locus).to_dict()

table = {'paper':   ['Sassetti 2003 in vivo','Lamichhane 2003',
                    'Sassetti 2003 in vitro','Rubin 2012'
                    ]
        }
df_table = pandas.DataFrame(table, columns = ['paper'])

def rename(name):
    if name == "E":
        name="Essential"
    elif name == "NE":
        name="Non-essential"
    else:
        name="-"
    return name

for gen in ernesto_FT_genes:
    data = [rename(sassetti2003Invivo_dic.get(gen,"-")), 
            rename(Lamichhane2003_dic.get(gen,"-")), 
            rename(sassetti2003Invitro_dic.get(gen,"-")), 
            rename(Rubin2012_dic.get(gen,"-"))]
    df_table[gen] = data
df_table= df_table.T
print(df_table)
df_table.to_csv('ogee_data.csv')