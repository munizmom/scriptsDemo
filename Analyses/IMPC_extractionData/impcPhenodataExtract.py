#!\usr\env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
#import os
import re
import datetime


'''
###################################################################
    Parsing script to extract relevant information from the 
    genotype-phenotype association data from the IMPC webpage
    in a tab based output table.

	Date: 25/03/19.

	Author: Mar Muñiz Moreno. PhD student Herault team @ IGBMC.

    Programming Language: python 2.7
###################################################################

README: This program needs two input arguments:
1. The formatted IMPC genotype-phenotype input data.
2. The name of the file where the results from this program will be written in. 
   If this file does not exist it will be created. If its exist, it will be append.

From the formatted IMPC genotype-phenotype input data (generation instruction over 
the script 
we will extract an output 
file formatted as a tab separated table compatible with downstream proccessing in 
excel and R. 
We will extract the following information:

a) markerSymbol: gene symbol.
b) numFound: number of phenotypes found associated with the gene.
c) mp_term_name: phenotypes terms.
d) top_level_mp_term_name: top level phenotype name, where the mp_terms are grouped inside.
e) intermediate_mp_term_name: interm level phenotype name, where the mp_terms are inside.
f) allele_symbol: mouse line name and allele.
g) parameter_name: the test or "parameter name" used where the phenotype was identified.
h) statistical_method: statistical test performed on the data collecter from the parameter.
i) percentage_change: % change between wt and allele mice observed. It can be + or - 
j) p_value: pvalue assigned to that change.
k) effect_size: effect size.
l) life_stage_name: related to the age of the mice analysed.
m) sex: sex of the mice where the phenotype was identified.
n) zygosity: zygosity state of the mice where the phenotype was identified.

How to run it:
Have activated python 2.7 in the environment.
As an example run the program using as input file Akt2_ex.txt by typing:
python scripts/impcPhenodataExtract.py trialData/Akt2_ex.txt trialData/output/trial.txt
python scripts/impcPhenodataExtract.py input/hsa21IMPC_formatted.txt resultsOutput/hsa21IMPC_dataExtracted.txt

'''


if __name__ == '__main__':
    input_filename = sys.argv[1]
    tab_name= sys.argv[2]
    try:
        #with open('/Users/marmmoreno/Documents/DS_models_book/ongoing/results/IMPC/dataRetrieval/hsa21IMPC_formatted.txt','r') as filename:
        with open(input_filename,'r') as filename:
            startTime = datetime.datetime.now()
            print("")
            print("") 
            print("          Procesing file ",filename)
            print("")
            print("          The parsing of the IMPC input began on ", str(startTime))
            print("")
            markerSymbol = ''
            numFound = ''
            mp_term_name = []
            top_level_mp_term_name = []
            intermediate_mp_term_name = []
            allele_symbol = []
            parameter_name = []
            statistical_method = []
            percentage_change = []
            p_value = []
            effect_size = []
            life_stage_name = []
            sex = []
            zygosity =[]


            with open(tab_name, 'a+') as tab_file: #add the header in output table
                header= "'{0}'\t'{1}'\t'{2}'\t'{3}'\t'{4}'\t'{5}'\t'{6}'\t'{7}'\t'{8}'\t'{9}'\t'{10}'\t'{11}'\t'{12}'\t'{13}'".format("Gene symbol", "Phenotypes Nb Found", "mp term name", "top level mp term name", "intermediate mp term name", "allele symbol", "sex", "zygosity", "parameter name", "statistical method", "percentage change", "p value", "effect size", "life stage name")
                tab_file.write(header + "\n")

            for line in filename:
                line = line.strip()
                if re.search('^}}', line) is not None :

                    with open(tab_name, 'a') as tab_file:
                        for x, y, z, v, k, oo, p, q, r, s, tt, u in map(None, mp_term_name, top_level_mp_term_name, intermediate_mp_term_name, allele_symbol,sex, zygosity, parameter_name, statistical_method, percentage_change, p_value, effect_size, life_stage_name):
                            sumUpTabInfo = "'{0}'\t'{1}'\t'{2}'\t'{3}'\t'{4}'\t'{5}'\t'{6}'\t'{7}'\t'{8}'\t'{9}'\t'{10}'\t'{11}'\t'{12}'\t'{13}'".format(markerSymbol, numFound, x, y, z, v, k, oo, p, q, r, s, tt, u)
                            tab_file.write(sumUpTabInfo + "\n")
                            #print (sumUpTabInfo)

                        markerSymbol = ''
                        numFound = ''
                        mp_term_name = []
                        top_level_mp_term_name = []
                        intermediate_mp_term_name = []
                        allele_symbol = []
                        parameter_name = []
                        statistical_method = []
                        percentage_change = []
                        p_value = []
                        effect_size = []
                        life_stage_name = []
                        sex = []
                        zygosity =[]


                else:

                    if re.search(r"\"q\":\"marker_symbol:" ,line):  #Search for the gene name
                        geneName = re.search(r"^\"q\":\"marker_symbol:(\w+)\"\,",line)
                        if geneName is not None :
                            markerSymbol = geneName.group(1) # Take gene_symbol name

                    elif re.search(r"^\"response\":\{\"numFound" ,line):
                        termsNb = re.search(r"^\"response\":\{\"numFound\":(\d+)\,.+" ,line)
                        if termsNb is not None :
                            numFound = termsNb.group(1)

                    elif re.search(r"^\"mp_term_name\":", line):
                        termMatch = re.search(r"^\"mp_term_name\":\"(.+)\"\,", line)
                        if termMatch is not None :
                            mp_term_name.append(termMatch.group(1))

                    # IMPORTANT NOTE1: for both top_level_mp_term_name & intermediate_mp_term_name 
                    # they can have so many terms that the IMPC decided to put it in 2 split lines
                    # or more. To know if is finished we need to find "]". That is why I read another
                    # line and try to find the end.
                    
                    # Note1 #####
                    elif re.search(r"^\"top_level_mp_term_name\":", line):
                        
                        while True:
                            if re.search(r"\]", line) is None:
                                line2=next(filename, '').strip()  
                                line = line + line2
                            else :
                                break
                            
                        termMatch = re.search(r"^\"top_level_mp_term_name\":\[\"(.+)\"\]\,", line)
                        if termMatch is not None:
                            top_level_mp_term_name.append(termMatch.group(1))
                            line=next(filename, '').strip()
                            if re.search(r"^\"intermediate_mp_term_id\":", line):
                                
                                while True:
                                    if re.search(r"\]", line) is None:
                                        line2=next(filename, '').strip()  
                                        line = line + line2
                                    else:
                                        break
                                line=next(filename, '').strip()

                                if re.search(r"^\"intermediate_mp_term_name\":", line):
                                    
                                    while True:
                                        if re.search(r"\]", line) is None:
                                            line2=next(filename, '').strip()  
                                            line = line + line2
                                        else:
                                            interTermMatch = re.search(r"^\"intermediate_mp_term_name\":\[\"(.+)\]\,", line)
                                            if interTermMatch is not None:
                                                intermediate_mp_term_name.append(interTermMatch.group(1))
                                            break
                                else:
                                    intermediate_mp_term_name.append('None')
                            
                            else:
                                intermediate_mp_term_name.append('None')

                            
                    #    interTermMatch = re.search(r"^\"intermediate_mp_term_name\":\[\"(.+)\]\,", line)
                    #    if interTermMatch is not None:
                    #        intermediate_mp_term_name.append(interTermMatch.group(1))
                            

                    ######
                    elif re.search(r"^\"allele_symbol\":", line):
                        alleleMatch = re.search(r"^\"allele_symbol\":\"(.+)\"\,", line)
                        if alleleMatch is not None :
                            allele_symbol.append(str(alleleMatch.group(1)))
     
                    elif re.search(r"^\"sex\":", line):
                        sexMatch = re.search(r"^\"sex\":\"(\w+)\"\,", line)
                        if sexMatch is not None :
                            sex.append(sexMatch.group(1))
          
                    elif re.search(r"^\"zygosity\":", line):
                        zygosMatch = re.search(r"^\"zygosity\":\"(\w+)\"\,", line)
                        if zygosMatch is not None :
                            zygosity.append(zygosMatch.group(1))


                    elif re.search(r"^\"parameter_name\":", line):
                        line=re.sub("\s\s+" , "", line) #remove the extra spaces found at the end of some lines inside the strings
                        paramMatch = re.search(r"^\"parameter_name\":\"(.+)\"\,", line)
                        if paramMatch is not None :
                            parameter_name.append(paramMatch.group(1))

                    elif re.search(r"^\"life_stage_name\":", line):
                        stageMatch = re.search(r"^\"life_stage_name\":\"(.+)\"\,", line)
                        if stageMatch is not None :
                            life_stage_name.append(stageMatch.group(1))

                    # NOTE IMPORTANT2: instead of being consistent with the formating,
                    # some entries may not have the 3 following lines after statistical_method:
                    # percentage_change, p_value and effect_size. 
                    # So we need to specify if they are not there, to make the values 'None'
                    # If not the next positive value will be wrongly stored here.
                    #So I need to make this recursive script below to account for that.

                    # Note2 #####
                    elif re.search(r"^\"statistical_method\":", line):                        
                        statsMatch = re.search(r"^\"statistical_method\":\"(.+)\"\,", line)
                        if statsMatch is not None :
                            statistical_method.append(statsMatch.group(1))
                        
                        line=next(filename, '').strip()  
                        
                        if re.search(r"^\"percentage_change\":", line):
                            changeMatch = re.search(r"^\"percentage_change\":\"(.+).\"\,", line)
                            if changeMatch is not None :
                                percentage_change.append(changeMatch.group(1))
                            else :
                                percentage_change.append('None')
                        
                            line=next(filename, '').strip()  
                            

                            if re.search(r"^\"p_value\":", line):
                                pvalMatch = re.search(r"^\"p_value\":(.+)\,", line)
                                if pvalMatch is not None :
                                    p_value.append(pvalMatch.group(1))
                                else :
                                    p_value.append('None')
                            
                                line=next(filename, '').strip()  

                                if re.search(r"^\"effect_size\":", line):
                                    eSizeMatch = re.search(r"^\"effect_size\":(.+)\,", line)
                                    if eSizeMatch is not None :
                                        effect_size.append(eSizeMatch.group(1))
                                    else :
                                        effect_size.append('None')
                                else :
                                    effect_size.append('None')

                            else :
                                p_value.append('None')
 
                        else :
                            percentage_change.append('None')
                            p_value.append('None')
                            effect_size.append('None')

                    ######
                #print(intermediate_mp_term_name)

        #formating the output in stout
        endTime = datetime.datetime.now() #printing the time of end
        print("          The parsing of the IMPC input finished on ", str(endTime), "took to run: ",str(endTime - startTime)) 
        print("") 
        print("          Script created by Mar Muniz Moreno. PhD student on 25/03/19")
        print("") 
        print("          --> Herault Lab @IGBMC <-- ") 


    #Raise exceptions in stout for the most common errors

    except IOError as e:
        print("          File reading error {0}: {1}".format(e.errno, e.strerror),file=sys.stderr)
    except:
        print("          Unexpected error: ", sys.exc_info()[0],file=sys.stderr)
        raise

else:
	raise AssertionError("          Please introduce at least one file with the correct IMPC genotype-phenotype data format.")
