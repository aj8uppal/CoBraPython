import os
import shutil

l0ed_list = [2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0]
fld_list = [4., 4.25, 4.5, 4.75, 5., 5.5, 6.]

for l0ed in l0ed_list:
    l1ed = l0ed
    for fld in fld_list:
        runname = 'l0ed_{}_l1ed_{}_fld_{}_highpower'.format(l0ed, l1ed, fld)
        os.system('python ./updateSizes.py -l0ed {} -l1ed {} -fld {}'.format(l0ed, l1ed, fld))
        os.system('python ./ITkBarrelManifold.py {}'.format(runname))

l0ed_list = [2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0]
fld_list = [4., 4.25, 4.5, 4.75, 5., 5.5, 6.]

for l0ed in l0ed_list:
    l1ed = l0ed
    for fld in fld_list:
        runname = 'l0ed_{}_l1ed_{}_fld_{}_highpower_2T'.format(l0ed, l1ed, fld)
        os.system('python ./updateSizes.py -l0ed {} -l1ed {} -fld {} -l0ec 4 -l1ec 4'.format(l0ed, l1ed, fld))
        os.system('python ./ITkBarrelManifold.py {}'.format(runname))

        
