import sys
import shutil
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument("-CRcl", dest="CRcaplength", type=float)
    parser.add_argument("-CRcd", dest="CRcapdiam", type=float)
    parser.add_argument("-CRec", dest="CRevpcorr", type=float)
    parser.add_argument("-CRvd", dest="CRevpdiam", type=float)
    parser.add_argument("-CR0ep", dest="CR0evppower", type=float)
    parser.add_argument("-CR1ep", dest="CR1evppower", type=float)

    parser.add_argument("-LRcl", dest="LRcaplength", type=float)
    parser.add_argument("-LRcd", dest="LRcapdiam", type=float)
    parser.add_argument("-LRec", dest="LRevpcorr", type=float)
    parser.add_argument("-LRvd", dest="LRevpdiam", type=float)
    parser.add_argument("-LRep", dest="LRevppower", type=float)

    parser.add_argument("-IRcl", dest="IRcaplength", type=float)
    parser.add_argument("-IRcd", dest="IRcapdiam", type=float)
    parser.add_argument("-IRec", dest="IRevpcorr", type=float)
    parser.add_argument("-IRvd", dest="IRevpdiam", type=float)
    parser.add_argument("-IRep", dest="IRevppower", type=float)

    parser.add_argument("-fld", dest="flexdiam", type=float)
    parser.add_argument("-exd", dest="exhdiam", type=float)
    
    args = parser.parse_args()

    lines = {}
    
    lines['CRcaplength'] = [18]
    lines['CRcapdiam'] = [10]
    lines['CR0evppower'] = [123,127]
    lines['CR1evppower'] = [125]
    lines['CRevpcorr'] = [139,141,143]
    lines['CRevpdiam'] = [11,12,13,14,15]
    
    lines['LRcaplength'] = [10]
    lines['LRcapdiam'] = [6]
    lines['LRevppower'] = [63]
    lines['LRevpcorr'] = [71]
    lines['LRevpdiam'] = [7]

    lines['IRcaplength'] = [10]
    lines['IRcapdiam'] = [6]
    lines['IRevppower'] = [63]
    lines['IRevpcorr'] = [71]
    lines['IRevpdiam'] = [7]

    lines['exhdiam'] = [5]

    lines['flexdiam'] = [5]
    
    coupledRing = ['../XML/QuarterShell/coupledRing_A_A.xml', '../XML/QuarterShell/coupledRing_A_B.xml', '../XML/QuarterShell/coupledRing_A_C.xml', '../XML/QuarterShell/coupledRing_A_D.xml',
                   '../XML/QuarterShell/coupledRing_B_A.xml', '../XML/QuarterShell/coupledRing_B_B.xml', '../XML/QuarterShell/coupledRing_B_C.xml', '../XML/QuarterShell/coupledRing_B_D.xml']
    coupledRing_vars = ['CRcaplength','CRcapdiam','CR0evppower','CR1evppower','CRevpcorr','CRevpdiam']

    l1Ring = ['../XML/QuarterShell/l1Ring_A_A.xml', '../XML/QuarterShell/l1Ring_A_B.xml', '../XML/QuarterShell/l1Ring_B_A.xml', '../XML/QuarterShell/l1Ring_B_B.xml']
    l1Ring_vars = ['LRcaplength','LRcapdiam','LRevppower','LRevpcorr','LRevpdiam']

    shortyRing = ['../XML/QuarterShell/shortyRing_A_A.xml', '../XML/QuarterShell/shortyRing_A_B.xml', '../XML/QuarterShell/shortyRing_B_A.xml', '../XML/QuarterShell/shortyRing_B_B.xml']
    shortyRing_vars = ['IRcaplength','IRcapdiam','IRevppower','IRevpcorr','IRevpdiam']

    flex = ['../XML/QuarterShell/flexlineA.xml', '../XML/QuarterShell/flexlineB.xml']
    flex_vars = ['flexdiam']

    exhaust = ['../XML/QuarterShell/exhaustA_A.xml', '../XML/QuarterShell/exhaustA_B.xml', '../XML/QuarterShell/exhaustB_A.xml', '../XML/QuarterShell/exhaustB_B.xml']
    exhaust_vars = ['exhdiam']

    sectors = [(coupledRing, coupledRing_vars), (l1Ring, l1Ring_vars), (shortyRing, shortyRing_vars), (exhaust, exhaust_vars), (flex, flex_vars)]
    
    for element, element_vars in sectors:
        for file in element:
            handle = open(file,'r')
            content = handle.readlines()
            for var in element_vars:
                if vars(args)[var] is not None:
                    for line in lines[var]:
                        content[line]='<value>{}</value>\n'.format(vars(args)[var])
            newfile = open('testEC.xml', 'w')
            [newfile.write(line) for line in content]
            newfile.close()
            handle.close()
            shutil.move('testEC.xml',  file)
        

