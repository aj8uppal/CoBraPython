import sys
import shutil
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument("-a", "--l0cl", dest="l0caplength", type=float)
    parser.add_argument("-b", "--l0cd", dest="l0capdiam", type=float)
    parser.add_argument("-c", "--l0ed", dest="l0exhdiam", type=float)
    parser.add_argument("-d", "--l0ec", dest="l0evpcorr", type=float)
    parser.add_argument("-e", "--l0ep", dest="l0evppower", type=float)
    
    parser.add_argument("-f", "--l1cl", dest="l1caplength", type=float)
    parser.add_argument("-g", "--l1cd", dest="l1capdiam", type=float)
    parser.add_argument("-h", "--l1ed", dest="l1exhdiam", type=float)
    parser.add_argument("-i", "--l1ec", dest="l1evpcorr", type=float)
    parser.add_argument("-j", "--l1ep", dest="l1evppower", type=float)
    
    parser.add_argument("-k", "--fld", dest="flexdiam", type=float)

    args = parser.parse_args()

    lines = {}
    
    lines['l0caplength'] = [20]
    lines['l0capdiam'] = [11]
    lines['l0exhdiam'] = [17]
    lines['l0evpcorr'] = [156,158,160]
    lines['l0evppower'] = [138,140,142]

    lines['l1caplength'] = [28]
    lines['l1capdiam'] = [15]
    lines['l1exhdiam'] = [25]
    lines['l1evpcorr'] = [198,200,202,204,206]
    lines['l1evppower'] = [224,226,228,230,232]

    lines['flexdiam'] = [5]
    
    staveL0 = ['../XML/HalfBarrel/staveL0_A.xml', '../XML/HalfBarrel/staveL0_B.xml']    
    staveL0_vars = ['l0caplength', 'l0capdiam', 'l0exhdiam', 'l0evpcorr', 'l0evppower']
    
    staveL1 = ['../XML/HalfBarrel/staveL1_A.xml', '../XML/HalfBarrel/staveL1_B.xml']
    staveL1_vars = ['l1caplength', 'l1capdiam', 'l1exhdiam', 'l1evpcorr', 'l1evppower']

    flex = ['../XML/HalfBarrel/flexlineA.xml', '../XML/HalfBarrel/flexlineB.xml']
    flex_vars = ['flexdiam']
    
    sectors = [(staveL0, staveL0_vars),(staveL1, staveL1_vars),(flex, flex_vars)]
    
    for element, element_vars in sectors:
        for file in element:
            handle = open(file,'r')
            content = handle.readlines()
            for var in element_vars:
                if vars(args)[var] is not None:
                    for line in lines[var]:
                        content[line]='<value>{}</value>\n'.format(vars(args)[var])
            newfile = open('test.xml', 'w')
            [newfile.write(line) for line in content]
            newfile.close()
            handle.close()
            shutil.move('test.xml',  file)
        

