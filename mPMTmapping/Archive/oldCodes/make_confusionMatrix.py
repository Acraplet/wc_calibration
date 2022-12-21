import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import read_data_results3 as rd
import matplotlib

plt.style.use(["science", "notebook", "grid"])
table = rd.read_data3("../table.txt")
ref = table[0]
print(ref)
compa = table[1]
print(compa)
LL = table[2]
LL_array = []
compa_array = []
ref_array = []


n = int(np.sqrt(len(ref)))
for j in range(n):
    buf_compa = []
    buf_LL = []
    buf_ref = []
    for i in range(int(np.sqrt(len(ref)))):
        print(ref[i+j*n], compa[i+j*n], LL[i+j*n])
        buf_ref.append(ref[i+j*n])
        buf_compa.append(compa[i+j*n])
        buf_LL.append(LL[i+j*n])
    print('\n')
    ref_array.append(buf_ref)
    compa_array.append(buf_compa)
    LL_array.append(buf_LL)

print(LL)

plt.figure(figsize = (10,10))
columns = []
for i in range(len(buf_compa)):
    columns.append("%i"%buf_compa[i])
#plt.xticks(ticks=np.arange(len(buf_ref)),labels=buf_ref,rotation=90)

hm = sns.heatmap(LL_array, cmap='viridis', norm=matplotlib.colors.LogNorm(), yticklabels=columns,xticklabels=columns, annot=True,  fmt=".2e")
#plt.colorbar(hm)
#plt.set_xticks(range(len(columns)),columns)
#plt.yticks(range(len(columns)),columns)
plt.title("Chi2 to the reference map")
plt.ylabel('Test map')
plt.xlabel('Reference map')
plt.show()
