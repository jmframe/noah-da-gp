#!/gpfsm/dulocal/sles11/other/SLES11.3/miniconda3/2019.03_py3.7/2019-05-15/bin/python
import numpy as np
import sys
#Compare recently compiled model code with the output of a tested Noah version.

with open('output_compare.out','r') as compareFile:
    compareModel = np.genfromtxt(compareFile)
with open('output.out','r') as newFile:
    newModel = np.genfromtxt(newFile)

if compareModel.shape != newModel.shape:
    print('The recent model run has a different output than the original')
    sys.exit()

for i in range(0,compareModel.shape[0]):
    for j in range(0,compareModel.shape[1]):
        if compareModel[i,j] != newModel[i,j]:
            print('The new model output does not match the original')
            print('at i = ', i, 'and j = ', j)
            print('original model output:')
            print(compareModel[i,:])
            print('new model output:')
            print(newModel[i,:])
            sys.exit()
print('The model output matches the original. Should be good to go')
