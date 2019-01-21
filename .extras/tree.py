import Imagetree as it
import numpy as np

mask = [[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,1,1,0],[0,0,0,0,0,1,1,0],[0,0,0,0,0,0,0,0]]
mask = np.array(mask)
print(mask)

for i in range(0,mask.shape[0]):
    for j in range(0,mask.shape[1]):
        print(mask[i][j])
    print('')

print(mask[0:mask.shape[0]/2,0:mask.shape[1]/2])
sum_val = np.empty(4)
#NE /SE/SW/NW

sum_val[0] = np.sum(mask[0:mask.shape[0]/2,mask.shape[1]/2:mask.shape[1]])
sum_val[1] = np.sum(mask[mask.shape[0]/2:mask.shape[0],mask.shape[1]/2:mask.shape[1]])
sum_val[2] = np.sum(mask[mask.shape[0]/2:mask.shape[0],0:mask.shape[1]/2])
sum_val[3] = np.sum(mask[0:mask.shape[0]/2,0:mask.shape[1]/2])
print(sum_val[0])
print(sum_val[1])
print(sum_val[2])
print(sum_val[3])

itrees =  it.Imagetree
for i in sum_val:
    if i == 0:



        