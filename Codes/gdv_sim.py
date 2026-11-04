import os
import sys
import os.path
import numpy as np
import math
import pandas as pd

selected = [1,2,3,5,6,7,8,9,10,11,12]

oi= [1.,2.,2.,3.,4.,3.,3.,4.,3.,4.,4.]
wi = [1.-(math.log(oi[i])/math.log(11.)) for i in range(11)]
swi = sum(wi)

def GDS_Distance(sig1, sig2):
    sdi = 0.
    for i in range(11):
        sdi += wi[i]*abs(math.log(sig1[i]+1.) - math.log(sig2[i]+1.))/(math.log(max(sig1[i], sig2[i])+2.))
    sdi /= swi
    return sdi

fname = "./Data/homo_sapiens_ppi_orca_counts.csv"
orca_counts_df = pd.read_csv(fname,index_col=0)
names = list(df.index)
GDVs = orca_counts_df.values

n = len(names)
print "Loaded %i_nodes x %i_orbits signatures"%(n, len(GDVs[0]))

full = (n*(n-1.))/2
perc = int(full/100.)
cnt = 0
DISTs = np.zeros((n,n))
for i in range(n):
	for j in range(i+1, n):
		cnt += 1
		if cnt%perc == 0:
			print "%i percent done"%(int(100.*cnt/float(full)))
		dij = GDS_Distance(GDVs[i], GDVs[j])
		DISTs[i][j] = 1-dij
		DISTs[j][i] = 1-dij

gdvd_df = pd.DataFrame(DISTs,index=names,columns=names)
gdvd_df.to_csv('./Data/gdvs_homo_sapiens_ppi.csv',header=True,index=True)
