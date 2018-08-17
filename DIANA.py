from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import copy
#

import sys
from Bio import SeqIO
from Bio import*
import os
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy as np
from Bio import AlignIO
from plotly.figure_factory._distplot import scipy
#from scipy.cluster import hierarchy
#from scipy.spatial.distance import pdist
import plotly.plotly
import plotly.figure_factory as ff
import math
np.set_printoptions(precision=5, suppress=True)

in_file = "dataset.txt"
seq_record =  SeqIO.parse(in_file , "fasta")
list_seq_record = list(seq_record)

#print(list_seq_record)
names=[]
for i in list_seq_record:
    names.append(i.id)


maxlen = max(len(record.seq) for record in list_seq_record)
#print(maxlen)
for record in list_seq_record:
    if len(record.seq) != maxlen :
        seq = str(record.seq).ljust(maxlen, '-')
        record.seq = Seq.Seq(seq)
assert all(len(record.seq)==maxlen for record in list_seq_record)

#out_file = '{}Padded.fasta'.format(os.path.splitext(in_file)[0])

with open("out_file", 'w') as fp:
    SeqIO.write(list_seq_record, fp, 'fasta')

alignment = AlignIO.read(open("out_file",'r'), 'fasta')

measure = DistanceCalculator('identity')
matrix = measure.get_distance(alignment)
#print matrix

pr = np.array(matrix)

#pr = np.array([[0.0,2.0,6.0,10.0,9.0],[2.0,0.0,5.0,9.0,8.0],[6.0,5.0,0.0,4.0,5.0],[10.0,9.0,4.0,0.0,3.0],[9.0,8.0,5.0,3.0,0.0]])

print (pr)
n = pr.shape[0]
print "Number of data points are " + str(n)

clus_matrix ={'LC':[],'RC':[],'dist':[],'size':[] }

current_cluster = []
dia_compare = [[],[]]

for steps in range(1,n):
    print "step number : " + str(steps)

    if(steps == 1 ):
        for i in range(0,n):
            current_cluster.append(i)
        #current_cluster=[0,1,2,3,4]

    print "dynamic current_cluster is" + str(current_cluster)
    dia = 0


    highest_mean = 0
    highest_mean_index =0
    for i in range(0,n):

        current_sum = 0

        for j in range(0,n):
            if(i in current_cluster and j in current_cluster):

                #below small 2 -line logic to handle case when all gives average zero.
                #So in this case highest_mean_index will point to 0,not to the first element(or any other element) of current_cluster
                #So in splinter group it will add element 0(bco we add element with highest_mean_index) in splinter
                #instead of adding first element(or any other element) of current cluster
                if(i == current_cluster[0]):
                    highest_mean_index = i

                if(pr[i][j] > dia):
                    dia = pr[i][j]
                if(j!=i):
                    current_sum =current_sum + pr[i][j]
        current_mean = current_sum / (len(current_cluster)-1)
        if(current_mean > highest_mean):
            highest_mean = current_mean
            highest_mean_index = i
    #print "Highest avg in this cluster is  " + str(highest_mean) + " at index " + str(highest_mean_index)

    splinter_group = []
    splinter_group.append(highest_mean_index)
    print splinter_group
    print str(steps) + "th step " + " splinter group contain " + str(splinter_group)

    non_splinter_group = []
    non_splinter_group = list(current_cluster)
    non_splinter_group.remove(highest_mean_index)
    print str(steps) + "th step " + " non - splinter group contain " + str(non_splinter_group)

    all_dif_neg_flag = 0
    while(all_dif_neg_flag == 0):
        highest_dif = -1000 #change it or logic after writing complete code
        for i in range(0, n):

            sum_non_splinter = 0
            sum_splinter = 0

            for j in range(0, n):
                if (i in non_splinter_group and j in non_splinter_group ):
                    if (j != i):
                        sum_non_splinter = sum_non_splinter + pr[i][j]

                if (i in non_splinter_group and j in splinter_group):
                    if (j != i):
                        sum_splinter =sum_splinter + pr[i][j]
            if (i in non_splinter_group): # we want to calculate difference for each non-splinter element
                #print "sum_non_splinter is " + str(sum_non_splinter)
                #print "sum_splinter is " + str(sum_splinter)


                if(sum_non_splinter != 0): #Handling singleton Case
                    avg_non_spliter = sum_non_splinter / (len(non_splinter_group) - 1)
                else:
                    avg_non_spliter = 0
                avg_spliter = sum_splinter / (len(splinter_group))

                avg_dif = avg_non_spliter - avg_spliter
                #print "avg difference for this element is  " + str(avg_dif)



                if (avg_dif > highest_dif):
                    highest_dif = avg_dif
                    highest_dif_index = i
        #end of for loop

        #print "element of non-splinter with highest-diff is at index " +str(highest_dif_index)


        if(highest_dif <= 0):
            splinter_group.sort()
            non_splinter_group.sort()
            clus_matrix['LC'].append(splinter_group)
            clus_matrix['RC'].append(non_splinter_group)
            clus_matrix['size'].append(len(splinter_group)+len(non_splinter_group))
            clus_matrix['dist'].append(dia)

            print "clus_matrix is : " + str(clus_matrix)

            print"Clustering of this step "+ str(steps)+" completed"
            break

        splinter_group.append(highest_dif_index)
        #print "Splinter group is " + str(splinter_group)

        non_splinter_group.remove(highest_dif_index)
        #print "non - Splinter group is " + str(non_splinter_group)

    #end of while loop

    #Calculate and compare diameter of splinter and non-splinter group and put the highest diameter group
    #into current_clus
    #but we need to compare with previous formed clusters also
    dia_splinter = 0
    dia_non_splinter = 0
    for i in range(0,n):

        for j in range(0,n):
            if(i in splinter_group and j in splinter_group):

                if(pr[i][j] > dia_splinter):
                    dia_splinter = pr[i][j]

    for i in range(0,n):

        for j in range(0,n):
            if(i in non_splinter_group and j in non_splinter_group):

                if(pr[i][j] > dia_non_splinter):
                    dia_non_splinter = pr[i][j]



    print "old  " + str(dia_compare)



    dia_compare[0].append(dia_non_splinter)
    dia_compare[1].append(non_splinter_group)

    dia_compare[0].append(dia_splinter)
    dia_compare[1].append(splinter_group)

    # When we have same 0 length sequences [28] and [119,120] bcoz of their inherent order [28] was coming before
    #[119,120] .So according to logic [28] will be retrieved instead of [119,128].It will affect logic.
    #So we sort everytime.Due to sorting inherent order of sequence merged is changed.
    # This is the reason that I am getting sorted dendogram.It is also correct.
    dia_compare[0].sort(reverse=True)
    dia_compare[1].sort(key=len,reverse=True)

    print "new " + str(dia_compare)

    max_dia = max(dia_compare[0])
    print max_dia

    max_dia_index = dia_compare[0].index(max_dia)

    current_cluster = dia_compare[1][max_dia_index]
    print current_cluster

    del dia_compare[0][max_dia_index]
    del dia_compare[1][max_dia_index]

    print dia_compare





print "Clustering Completed.Final clus_matrix is "
print clus_matrix

index_array = []

Z_list = [] #[[],[],[],[]] #it should be dynamic
#print Z_list[1]
for i in range(0,n):
    index_array.append([i])

print "index array initially is " +str(index_array)


for i in range(n-2,-1,-1):
    print clus_matrix['LC'][i]
    print clus_matrix['RC'][i]
    Z_list_ka_index = (n - 2) - i
    index_array_ka_index = n + Z_list_ka_index #(2*n -2) - i
    print "Searching "+ str(i) + " entry of clus_matrix"
    print "This step will fill "+ str(index_array_ka_index) + "th index of index_array and fill "+ str(Z_list_ka_index) + "th list of Z_array"
    temp1_list = []
    temp2_list = []
    Z_list.append([])
    for j in range((index_array_ka_index -1),-1,-1 ):
        if(clus_matrix['LC'][i] == index_array[j]):

            Z_list[Z_list_ka_index].append(j)
            temp1_list = list(clus_matrix['LC'][i])
            #temp1_list.sort()

        if(clus_matrix['RC'][i] == index_array[j]):

            Z_list[Z_list_ka_index].append(j)
            temp2_list = list(clus_matrix['RC'][i])
            #temp2_list.sort()
    merged_list = temp1_list + temp2_list
    merged_list.sort()
    index_array.append(merged_list)
    print "now index array is  " + str(index_array)
    Z_list[Z_list_ka_index].append(clus_matrix['dist'][i])
    Z_list[Z_list_ka_index].append(clus_matrix['size'][i])
#Z_list[2] = list(clus_matrix['dist'])
#Z_list[3] = list(clus_matrix['size'])
print "index array is " +str(index_array)
print "Z_list is " +str(Z_list)

Z=np.array(Z_list)
print "Z array is " +str(Z)

plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()











