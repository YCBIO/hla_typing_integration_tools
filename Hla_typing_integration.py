#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys, os, argparse, re, math, glob, subprocess, copy
import getopt

def USAGE(script):
    helpMess = '''
     ProgramName:\t%s
         Version:\tV1.0
         Contact:\tyijian@yucebio.com
    Program Date:\t2019.09.11
          Modify:\t-
     Description:\tThis program is used to integration HLA typing result of 
                   optitype and polysolver.
           Usage:\tpython3 %s -o <optityperesult>  -p <polyresult>  -r <result> 
         Options:
            -h               Show this help message.
            -o --optiresult  optitype typing result
            --on             optitype typing result of normal sample 
            --ot             optitype typing result of tumor sample
            --pn             polysolver typing result of normal sample 
            --pt             polysolver typing result  of tumor sample 
            -p --polyresult  polysolver typing result
            -r --result      output result file(required)
    ''' %(script,script)
    return helpMess


def read_optitype(optires):
    hladict = {}
    linenum = 1
    for line in open(optires,"r"):
        if linenum == 2:
            linelist = line.split('\t')
            alleles = []
            for  a in [linelist[1],linelist[2]]:
                if a != "":
                    alleles.append(a)
            hladict["HLA-A"] = alleles
            alleles = []
            for  a in [linelist[3],linelist[4]]:
                if a != "":
                    alleles.append(a)
            hladict["HLA-B"] = alleles
            alleles = []
            for  a in [linelist[5],linelist[6]]:
                if a != "":
                    alleles.append(a)
            hladict["HLA-C"] = alleles
            #hladict = { "HLA-A":[linelist[1],linelist[2]],"HLA-B":[linelist[3],linelist[4]],"HLA-C":[linelist[5],linelist[6]]}
            #print(hladict)
            return hladict
        linenum += 1


def read_polysolver(polyres):
    hladict = {}
    genelist = ["HLA-A","HLA-B","HLA-C"]
    linenum = 0
    for line in open(polyres,"r"):
        line = line.strip('\n')
        linelist = line.split('\t')
        #print(linelist)
        if linelist[1]:
            allele1list = linelist[1].split('_')
            allele1 = allele1list[1].upper()+"*"+allele1list[2]+":"+allele1list[3]
        else:
            allele1 = ""

        if linelist[2]:
            allele2list = linelist[2].split('_')
            allele2 = allele2list[1].upper()+"*"+allele2list[2]+":"+allele2list[3]
        else:
            allele2 = ""
        #hlalist.append([allele1,allele2])
        alleles = []
        for a in [allele1,allele2]:
            if a != "":
                alleles.append(a)
        hladict[genelist[linenum]] = alleles
        linenum += 1
    #print(hladict)
    return hladict

def integration_hlares(opti_n,opti_t,poly_n,poly_t,result):
    outfile = open(result,"w")
    
    hla_opti_n = copy.deepcopy(read_optitype(opti_n))
    hla_opti_t = copy.deepcopy(read_optitype(opti_t))
    hla_poly_n = copy.deepcopy(read_polysolver(poly_n))
    hla_poly_t = copy.deepcopy(read_polysolver(poly_t))
    hla_opti_n_string = creat_hla_res_string(hla_opti_n)
    hla_opti_t_string = creat_hla_res_string(hla_opti_t)
    hla_poly_n_string = creat_hla_res_string(hla_poly_n)
    hla_poly_t_string = creat_hla_res_string(hla_poly_t)
    integration = ""
    for g in ["HLA-A","HLA-B","HLA-C"]:
        hla_opti_n_g = hla_opti_n[g]
        hla_opti_t_g = hla_opti_t[g]
        hla_poly_n_g = hla_poly_n[g]
        hla_poly_t_g = hla_poly_t[g]
        if set(hla_opti_n_g) == set(hla_opti_t_g):
            [allele1,allele2] = sorted(hla_opti_n_g)
        else:
            [allele1,allele2] = sorted(count_score(hla_opti_n_g,hla_opti_t_g,hla_poly_n_g,hla_poly_t_g))
        integration += "HLA-%s,HLA-%s," % (allele1,allele2)
    integration = integration.strip(",")
    outfile.write("Integration\t%s\n" % (integration))
    outfile.write("Optitype_n\t%s\nOptitype_t\t%s\nPolysolver_n\t%s\nPolysolver_t\t%s\n" % (hla_opti_n_string,hla_opti_t_string,hla_poly_n_string,hla_poly_t_string))

def creat_hla_res_string(resDict):
    resString = ""
    for g in ["HLA-A","HLA-B","HLA-C"]:
        for h in sorted(resDict[g]):
            resString += "HLA-%s," % (h)
    resString = resString.strip(",")
    return resString

def count_score(opti_n,opti_t,poly_n,poly_t):
    Weights = {"opti_n":100,"opti_t":98,"poly_n":96,"poly_t":91}
    scoredict = {}
    for g in opti_n:
        key = g
        if g in scoredict.keys():
            key = g + "_pure"
        scoredict[key] =  Weights["opti_n"]
        if g in opti_t:
            scoredict[key] += Weights["opti_t"]
            opti_t.remove(g)
        if g in poly_n:
            scoredict[key] += Weights["poly_n"]
            poly_n.remove(g)
        if g in poly_t:
            scoredict[key] += Weights["poly_t"]
            poly_t.remove(g)
    for g in opti_t:
        key = g
        if g in scoredict.keys():
            key = g + "_pure"
        scoredict[key] =  Weights["opti_t"]
        if g in poly_n:
            scoredict[key] += Weights["poly_n"]
            poly_n.remove(g)
        if g in poly_t:
            scoredict[key] += Weights["poly_t"]
            poly_t.remove(g)
    for g in poly_n:
        key = g
        if g in scoredict.keys():
            key = g + "_pure"
        scoredict[key] =  Weights["poly_n"]
        if g in poly_t:
            scoredict[key] += Weights["poly_t"]
            poly_t.remove(g)
    for g in poly_t:
        key = g
        if g in scoredict.keys():
            key = g + "_pure"
        scoredict[key] =  Weights["poly_t"]
    scoredict_sorted_list = sorted(scoredict.items(), key=lambda item: item[1], reverse=True)
    allele1 = scoredict_sorted_list[0][0]
    allele2 = scoredict_sorted_list[1][0]
    allele2 = allele2.split("_pure")[0]
    return [allele1,allele2]

if __name__ == '__main__':
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"ho:p:r:",["optiresult=","polyresult=","result=","on=","ot=","pn=","pt="])
    except getopt.GetoptError:
        print(USAGE(sys.argv[0]))
        sys.exit(2)
    if len(opts) == 0:
        print(USAGE(sys.argv[0]))
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h','--help'):
            print(USAGE(sys.argv[0]))
            sys.exit()
        elif opt in ("-o", "--optitype"):
            optiresult = arg
            optiresult = os.path.abspath(optiresult)
            if not os.path.exists(optiresult):
                print ('ERROR: %s dose not exists!' % optiresult)
                sys.exit(1)
        elif opt in ("-p", "--polyresult"):
            polyresult = arg
            polyresult = os.path.abspath(polyresult)
            if not os.path.exists(polyresult):
                print ('ERROR: %s dose not exists!' % polyresult)
                sys.exit(1)
        elif opt in ( "--on"):
            optinormal = arg
            optinormal = os.path.abspath(optinormal)
            if not os.path.exists(optinormal):
                print ('ERROR: %s dose not exists!' % optinormal)
                sys.exit(1)
        elif opt in ( "--ot"):
            optitumor = arg
            optitumor = os.path.abspath(optitumor)
            if not os.path.exists(optitumor):
                print ('ERROR: %s dose not exists!' % optitumor)
                sys.exit(1)
        elif opt in ( "--pn"):
            polynormal = arg
            polynormal = os.path.abspath(polynormal)
            if not os.path.exists(polynormal):
                print ('ERROR: %s dose not exists!' % polynormal)
                sys.exit(1)
        elif opt in ( "--pt"):
            polytumor = arg
            polytumor = os.path.abspath(polytumor)
            if not os.path.exists(polytumor):
                print ('ERROR: %s dose not exists!' % polytumor)
                sys.exit(1)
        elif opt in ("-r", "--result"):
            result = arg
            result = os.path.abspath(result)
    integration_hlares(optinormal,optitumor,polynormal,polytumor,result)
