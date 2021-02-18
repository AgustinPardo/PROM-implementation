iEK1011_file=open("iEK1011_genes.txt","r")
iEK1011_genes=[i.rstrip() for i in iEK1011_file.readlines()]
iEK1011_genes_set= set(iEK1011_genes)

targets_file=open("targets.txt","r")
targets_genes=[i.rstrip() for i in targets_file.readlines()]
targets_genes_set = set(targets_genes)

targets_inter_iEK1011=targets_genes_set-iEK1011_genes_set

#Leo la red
red=open("ernesto_net.txt","r")
red_lines=red.readlines()

target_filter=open("target_filter.txt","w")
regulator_filter=open("regulator_filter.txt","w")
litevidence_filter=open("litevidence_filter.txt","w")
prob_prior_filter=open("prob_prior_filter.txt","w")

for line in red_lines:
    regulator=line.split(",")[0]
    target=line.split(",")[1].rstrip()
    if regulator == target:
        pass
    else:
        if target in targets_inter_iEK1011:
            pass
        else:
            target_filter.write(target+"\n")
            regulator_filter.write(regulator+"\n")
            litevidence_filter.write("0"+"\n")
            prob_prior_filter.write("0"+"\n")

iEK1011_file.close()
targets_file.close()
red.close()
target_filter.close()
regulator_filter.close()
litevidence_filter.close()
prob_prior_filter.close()