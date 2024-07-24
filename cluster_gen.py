import pandas as pd
from Bio import Align
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans



ic50=pd.read_csv("/Users/jasonhwang/Documents/heat_map/ic50_cov_Val.csv")
abs_name=ic50["Antibody  Name"]
vh=ic50["Heavy chain AA"]
vl=ic50["Light chain AA"]
gene_vl=ic50["Light chain V gene"]
gene_vl=ic50["Heavy chain V gene"]

sources=ic50["source"]
print(len(sources))
sorc=[]
for i in range(len(sources)):
    sorc.append("6")
for i in range(len(sources)):
    if "BA.1 convalescents" in sources[i]:
        sorc[i]="0"
    if "BA.2 convalescents" in sources[i]:
        sorc[i]="1"
    if "BA.5 convalescents" in sources[i]:
        sorc[i]="2"
    if "WT convalescents" in sources[i]:
        sorc[i]="3"
    if "SARS convalescents" in sources[i]:
        sorc[i]="4"
    if "WT vaccinees" in sources[i]:
        sorc[i]="5"



cluster_val=vh
clus_align=[]

clus_sim=[]
def seq_aligner(seq2):
    # Create a dictionary of sequences
    # Initialize the aligner
    aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
    aligner.mode = 'local'
    
    # Reference sequence
    seq1 = "EVQLVESGGGLIQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"
    
    # Align seq2 with the reference sequence seq1
    alignments = aligner.align(seq1, seq2)
    
    # Get the best alignment (first alignment in the list)
    best_alignment = alignments[0]
    # Store the sequence and its alignment
    clus_align.append(best_alignment[1])
    tempscore=aligner.score(seq1, seq2)
    clus_sim.append(tempscore)
    #print(tempscore)
    #for line in best_alignment:
        #print(len(line))
     #   vh_align.append(line)
        
for i in range(len(cluster_val)):
    seq_aligner(cluster_val[i])
cluster_convert=[]
cluster_convert_val=[]
for i in range(len(clus_align)):
    #print(len(clus_align[i]))
    cluster_convert.append(clus_align[i])
    value=0
    for j in range(len(clus_align[i])):
        #print(vh_align[i][j])
        x=ord(clus_align[i][j])
        #print(x)
        value+=x
    cluster_convert_val.append(value)
#x=cluster_convert_val
x=clus_sim
#plt.scatter(x, y)
#plt.show()

g=open("x.txt", "w")
for line in clus_align:
    g.write(f"{line} \n")

######
from sklearn.cluster import KMeans

name_ascii=[]
for i in range(len(abs_name)):
    count=0
    for j in range(len(abs_name[i])):
        z=ord(abs_name[i][j])
        count+=z
    name_ascii.append(count)
name_ascii_max=max(name_ascii)
for i in range(len(name_ascii)):
    name_ascii[i]=name_ascii[i]/name_ascii_max*10

name_ascii_max=max(name_ascii)
y=name_ascii
data = list(zip(x, y))
inertias = []

for i in range(1,11):
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(data)
    inertias.append(kmeans.inertia_)

plt.plot(range(1,11), inertias, marker='o')
#plt.title('Elbow method')
#plt.xlabel('Number of clusters')
#plt.ylabel('Inertia')

kmeans = KMeans(n_clusters=6)
kmeans.fit(data)
plt.scatter(x, y, c=kmeans.labels_)
plt.show()

cluster_label=kmeans.labels_

BG=ic50["D614G"]
for i in range(len(BG)):
    if BG[i]==">10":
        BG[i]="10"
BA1=ic50["BA.1"]
for i in range(len(BA1)):
    if BA1[i]==">10":
        BA1[i]="10"
BA2=ic50["BA.2"]
for i in range(len(BA2)):
    if BA2[i]==">10":
        BA2[i]="10"
BA5=ic50["BA.5"]
for i in range(len(BA5)):
    if BA5[i]==">10":
        BA5[i]="10"
BA275=ic50["BA.2.75"]
for i in range(len(BA275)):
    if BA275[i]==">10":
        BA275[i]="10"
BA2752=ic50["BA.2.75.2"]
for i in range(len(BA2752)):
    if BA2752[i]==">10":
        BA2752[i]="10"
CA1=ic50["CA.1"]
for i in range(len(CA1)):
    if CA1[i]==">10":
        CA1[i]="10"
BQ11=ic50["BQ.1.1"]
for i in range(len(BQ11)):
    if BQ11[i]==">10":
        BQ11[i]="10"
BR2=ic50["BR.2"]
for i in range(len(BR2)):
    if BR2[i]==">10":
        BR2[i]="10"
BM=ic50["BM.1.1.1"]
for i in range(len(BM)):
    if BM[i]==">10":
        BM[i]="10"
XBB=ic50["XBB"]
for i in range(len(XBB)):
    if XBB[i]==">10":
        XBB[i]="10"
cluster_data=pd.DataFrame({
    "name":abs_name, 
    "name_ascii code":name_ascii, 
    "vh seq":vh, 
    "source clustering":sorc, 
    "cluster":cluster_label, 
    "D614G":BG, 
    "BA.1":BA1, 
    "BA.2":BA2, 
    "BA.5":BA5, 
    "BA.2.75":BA275, 
    "BA.2.75.2":BA2752, 
    "CA.1":CA1, 
    "BQ.1.1":BQ11, 
    "BR.2":BR2, 
    "BM.1.1.1":BM, 
    "XBB":XBB
    })

sorted_cluster=cluster_data.sort_values(by="cluster", ascending=True)
sorted_cluster.to_csv("ic50_clustered_data.csv")