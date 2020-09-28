import os
import shutil

exd_list = [4.0, 5.0, 6.0, 7.0]
fld_list = [4.0, 5.0, 6.0, 7.0]

for exd in exd_list:
    for fld in fld_list:
        runname = 'exd_{}_fld_{}_highpowerEC'.format(exh, fld)
        os.system('python ./updateSizesEC.py -exd {} -fld {}'.format(exd, fld))
        os.system('python ./ITkEndCapManifold.py {}'.format(runname))

exd_list = [4.0, 5.0, 6.0, 7.0]
fld_list = [4.0, 5.0, 6.0, 7.0]

for exd in exd_list:
    for fld in fld_list:
        runname = 'exd_{}_fld_{}_highpowerEC_2T'.format(exh, fld)
        os.system('python ./updateSizesEC.py -exd {} -fld {} -CRec 4 -LRec 4 -IRec 4'.format(exd, fld))
        os.system('python ./ITkEndCapManifold.py {}'.format(runname))
