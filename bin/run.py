import os
import shutil

l0cd_list = [0.85, 0.8, 0.75, 0.7, 0.65]
l1cd_list = [1.25, 1.2, 1.15, 1.10, 1.05]
l0cl_list = [2.7, 2.6, 2.5, 2.4, 2.3, 2.2]
l1cl_list = [2.7, 2.6, 2.5, 2.4, 2.3, 2.2]

for l0cd in l0cd_list:
    for l1cd in l1cd_list:
        for l0cl in l0cl_list:
            for l1cl in l1cl_list:
                runname = 'l0cd_{}_l1cd_{}_l0cl_{}_l1cl_{}'.format(l0cd, l1cd, l0cl, l1cl)
                os.system('python ./updateSizes.py --l0cd {} --l1cd {} --l0cl {} --l1cl {}'.format(l0cd,l1cd,l0cl,l1cl))
                os.system('python ./ITkBarrelManifold.py {}'.format(runname))
