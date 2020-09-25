import os
import shutil

l0ed_list = [4., 4.5, 5., 5.5, 6.]
fld_list = [4., 4.5, 5., 5.5, 6.]

for l0ed in l0ed_list:
    l1ed = l0ed
    for fld in fld_list:
        runname = 'l0ed_{}_l1ed_{}_fld_{}_highpower'.format(l0ed, l1ed, fld)
        os.system('python ./updateSizes.py -l0ed {} -l1ed {} -fld {}'.format(l0ed, l1ed, fld))
        os.system('python ./ITkBarrelManifold.py {}'.format(runname))

