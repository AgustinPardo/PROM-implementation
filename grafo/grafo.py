import networkx as nx
import matplotlib.pyplot as plt
import collections
import numpy as np

FT=open("/home/agustin/FBA_Tesis/trabajo/regulator_filter.txt","r")
target=open("/home/agustin/FBA_Tesis/trabajo/targets_filter.txt","r")

clean = lambda lista : [i.rstrip() for i in lista]

FT_lines=clean(FT.readlines())
target_lines=clean(target.readlines())

FT_set=set(FT_lines)
target_set=set(target_lines)

print("#FT:"+  str(len(FT_set)))
print("#target:"+  str(len(target_set)))
print("#FT & target:"+ str(len(FT_set & target_set)))

DG = nx.DiGraph()


#DG.add_edges_from([(3, 4), (4, 5)], color="red")
for index in range(len(FT_lines)):
    DG.add_weighted_edges_from([(FT_lines[index], target_lines[index], 1)])

options = {
    'with_labels' : True,
}

# nx.draw(DG, **options)
# plt.show()

alld=(nx.info(DG))
ind=(DG.in_degree())
outd=(DG.out_degree())

ind = [i for i in ind if i[1] >0]
outd = [i for i in outd  if i[1] > 0]

def degree_hist(degree_dist):
    degree_sequence = sorted([d for n, d in degree_dist], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())   
    return [deg, cnt]


fig, (ax1, ax2) = plt.subplots(2)
ax1.bar(degree_hist(outd)[0],degree_hist(outd)[1], width=0.80, color="b")

ax1.set_title("Out Degree Histogram")
ax1.set(xlabel="Degree",ylabel="count")
ax1.set_xticks(np.arange(0, max(degree_hist(outd)[0]), step=5))

ax2.bar(degree_hist(ind)[0],degree_hist(ind)[1], width=0.80, color="b")
ax2.set_title("In Degree Histogram")
ax2.set(xlabel="Degree",ylabel="count")
ax2.set_xticks(np.arange(0, max(degree_hist(ind)[0]), step=5))

fig.tight_layout()
plt.show()