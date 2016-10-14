# must common
import os
import sys
import itertools
import math

# for regular expressions
import re

# To manipulate spreadsheets of excel
import openpyxl
import xlsxwriter

# from Biopython
from Bio import GenBank
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Other Important Biopython modules
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Blast.Applications import NcbiblastnCommandline


# Functions not mine

############################################################

def Processing_files_from_folder(path):
    
    dirExists = os.path.isdir(path)
    if(dirExists != False):
        # como no se sabe el numero de genomas a ingresar
        list_files = []
        for(files) in os.walk(path):
        ## printing files to see 
            list_files = files
        
        list_genomes_names = list_files[2] 
    
        for geno in list_genomes_names:
            if(geno.startswith('._') or geno.startswith('.DS_') ):
                list_genomes_names.remove(geno)
        
        mssge = "Folder exists"
    else:
        list_genomes_names = []
        mssge = "Folder does not exist"

    return(list_genomes_names,mssge)

def Processing_folder_paths(path):
    
    i = 0
    folders = []
    for folder in os.walk(path):

        i += 1

        if(i == 1):
            mssg = "All folder paths were obtained"
            folders = folder[1]
            break
            
    if(len(folders) == 0):
        mssg = "There are not folder in this PATH"
        
    print mssg
    
    return folders

def Getting_genbank_files_names(list_genomes_names):
    
    names = []
    gb_files = []
    
    # Doing loop to process fasta files
    for f in list_genomes_names:
        
        if(f.endswith(".gb")):
        # Names for outpfiles
            base = os.path.splitext(f)[0]
            names.append(base)
        
            gb_files.append(f)
     
    return gb_files, names

def Parsing_fasta(fileName):
    
    handleMultifasta = open(fileName,"rU")
    
    fileFasta = SeqIO.parse(handleMultifasta,"fasta")
    
    return fileFasta

def fasta_reader(filename):
    
    from Bio.SeqIO.FastaIO import FastaIterator
    with open(filename) as handle:
    
        for record in FastaIterator(handle):
        
            yield record

def Parsing_genbank_with_SeqIO(fileName):
    
    handleMultifasta = open(fileName,"r")
    
    fileGenbank = SeqIO.parse(handleMultifasta,"genbank")
    
    return fileGenbank

def Parsing_gb_with_GenBank(fileName):
    
    from Bio import GenBank
    
    parser = GenBank.RecordParser()
    record = parser.parse(open(fileName))
    
    return record

# Making the openpyxl object to manipulate each .xlsx where the scores are
def Making_sheet_for_overwrite_files_active(path):
    
    # Creating workbook object of the .xlsx
    spSheetWrite = openpyxl.load_workbook(path)

    # Activating this workbook to be overwrited
    spSheetWriteActive = spSheetWrite.active
    
    return(spSheetWriteActive, spSheetWrite)

def Getting_files_path_by_extension(list_genomes_names,ext):
    
    names = []
    path_files = []
    
    # Doing loop to process fasta files
    for f in list_genomes_names:
        
        if(f.endswith(str(ext))):
        # Names for outpfiles
            base = os.path.splitext(f)[0]
            names.append(base)
        
            path_files.append(f)
     
    return path_files, names

## For mapping with blastn two seq
def magnitude(x):
    import math
    "Not mine but it is very simple"
    return int(math.floor(math.log10(x)))

##############################################################
### Functions de Ribosomal Operons Multiple Files To Excel => First Properties
##############################################################

## Auxiliar Function 0.0
def unique_list_function(list_sample):
    
    unique_list = [e for i, e in enumerate(list_sample) if list_sample.index(e) == i]
    
    return unique_list

# Auxiliar Function 0.1
def getting_operon_count_in_gb_order(uniqueNamesList,operon_count_names,operon_count_values):
    
    countOrderedAsUnique = []
    
    for uniqueName in uniqueNamesList:
        
        for name , value in itertools.izip(operon_count_names, operon_count_values):
            
            if(uniqueName == name):
                countOrderedAsUnique.append(value)
                
    return countOrderedAsUnique

## Function 0
def getting_count_each_operon_2_list(operonNamesTemp):
    
    uniqueNamesList = unique_list_function(operonNamesTemp)

    operon_count_dict = {rna:operonNamesTemp.count(rna) for rna in operonNamesTemp}
    operon_count_names = operon_count_dict.keys()
    operon_count_values = operon_count_dict.values()
    
    countOrderedAsUnique = getting_operon_count_in_gb_order(uniqueNamesList,operon_count_names,operon_count_values)

    return uniqueNamesList, countOrderedAsUnique

## Function 1 
def making_dict_by_geno(gb_generator_list):
    
    dict_geno_gb_properties = {}


    for geno in gb_generator_list:

        for f in geno:

            #print f.id
            nameT = f.name
            #namesGenos.append(nameT)

            #print nameT

            listTempGeno = []

            #creatign RNA list
            rRNAs=[]
            operonNames = []

            #other lists
            product_list=[]
            product_split_final_list = []
            #location_list=[]
            start_list=[]
            end_list=[]
            strand_list=[]
            desc_list=[]
            resta_list=[]

            for feature in f.features:

                if(feature.type == "rRNA"):

                    start = feature.location.start
                    end = feature.location.end
                    #location = str(start) + ":" + str(end)
                    resta = end - start

                    desc = feature.qualifiers['locus_tag'][0]            
                    product = feature.qualifiers['product'][0]

                    product_split = product.split(" ")
                    #print len(product_split)

                    if(len(product_split) >= 5):
                        product_split_final = product_split[4]
                    else:
                        product_split_final = "No_def"

                    ## appending to product_split_final
                    product_split_final_list.append(product_split_final)


                    seq = feature.extract(f.seq)
                    strand_dir = feature.location.strand


                    #lists
                    product_list.append(product)
                    #location_list.append(location)
                    operonNames.append(product_split_final)
                    start_list.append(start)
                    end_list.append(end)
                    strand_list.append(str(strand_dir))
                    desc_list.append(desc)
                    resta_list.append(resta)
                    
            count_operon_names, count_operon_list = getting_count_each_operon_2_list(operonNames)
            dict_geno_gb_properties[nameT] = (product_list, start_list, end_list, strand_list, resta_list, count_operon_names,
                                             count_operon_list)
            
    return(dict_geno_gb_properties)

## Function 1 Version 2
def making_dict_by_geno_from_multiGenbank(genbankParser):
    
    dict_geno_gb_properties = {}

    for f in genbankParser:
        
        nameT = f.name
            
        listTempGeno = []

        rRNAs=[]
        operonNames = []

        product_list=[]
        product_split_final_list = []
               
        start_list=[]
        end_list=[]
        strand_list=[]
        desc_list=[]
        resta_list=[]

        for feature in f.features:

                if(feature.type == "rRNA"):

                    start = feature.location.start
                    end = feature.location.end
                        
                    resta = end - start

                    desc = feature.qualifiers['locus_tag'][0]            
                    product = feature.qualifiers['product'][0]

                    product_split = product.split(" ")
                        

                    if(len(product_split) >= 5):
                        product_split_final = product_split[4]
                    else:
                        product_split_final = "No_def"

                    ## appending to product_split_final
                    product_split_final_list.append(product_split_final)


                    seq = feature.extract(f.seq)
                    strand_dir = feature.location.strand


                    #lists
                    product_list.append(product)
                    operonNames.append(product_split_final)
                    start_list.append(start)
                    end_list.append(end)
                    strand_list.append(str(strand_dir))
                    desc_list.append(desc)
                    resta_list.append(resta)

        count_operon_names, count_operon_list = getting_count_each_operon_2_list(operonNames)
        dict_geno_gb_properties[nameT] = (product_list, start_list, end_list, strand_list, resta_list, count_operon_names,
                                                 count_operon_list)

    return(dict_geno_gb_properties)

## Function 2
def getting_blastn_cline2_rRNA_vs_geno_getting_positions(pathMain,operon_list,fasta_files,uniqueGenoNames):
    
    dict_blast_split_lines = {}
    
    for operonfile, fastafile, uniqueGenoName in itertools.izip(operon_list,fasta_files,uniqueGenoNames):
        
        name = uniqueGenoName
        
        finalPath = pathMain + operonfile
        fastaPath = pathMain + fastafile
    
        # using NcbiblastnCommandline
        blastn_cline2 = NcbiblastnCommandline(query= finalPath, 
                                         subject = fastaPath, 
                                         outfmt = 6, max_hsps = 2)()[0]
        
        blastn_cline2_split_lines = blastn_cline2.splitlines( )
    
        dict_blast_split_lines[name] = blastn_cline2_split_lines
        
    return dict_blast_split_lines

## Function 3
def making_map1(blastn_cline2_split_lines):
    
    mapping_locations_geno2 = []
    list_start_map_geno2 = []
    list_end_map_geno2 = []

    for line in blastn_cline2_split_lines:
    
        location = []
    
        ident_seq =line.split("\t")[0]
    
        start = line.split("\t")[8]
        end = line.split("\t")[9]
    
        list_start_map_geno2.append(start)
        list_end_map_geno2.append(end)
    
        position = start+":"+end
    
        location.extend((ident_seq,position))
    
        mapping_locations_geno2.append(location)
        
    return(mapping_locations_geno2)

## Function 4
def using_making_map1_on_dict(dict_from_blastn_cline2):
    
    dict_using_making_map1 = {}
    
    for key in dict_from_blastn_cline2:
        
        blastn_splitlinesTemp = dict_from_blastn_cline2[key]
        
        making_map1_resultTemp = making_map1(blastn_splitlinesTemp)
        
        dict_using_making_map1[key] = making_map1_resultTemp 
        
    return dict_using_making_map1

## Function 5
def making_map2(mapping_locations_geno2):
    
    import collections
    
    #mapping_locations_geno2
    mapping_locations_geno2 = tuple(mapping_locations_geno2)

    data = collections.defaultdict(list)
    
    for k, v in mapping_locations_geno2:
    
        data[k].append(v)
    
    return data

## Function 6
def using_making_map2_on_dict(dictMap1):
    
    dict_using_making_map2 = {}
    
    for key in dictMap1:
        
        maplocationsMap1 = dictMap1[key]
        
        data_map2 = making_map2(maplocationsMap1)
        
        dict_using_making_map2[key] = data_map2
    
    return dict_using_making_map2

## Function 7
def gettting_position_with_less_o_magnitude(data):

    mapping_locations_geno3 = {}

    for key,value in data.items():
    
    #print key, value
    
        primero = value[0].split(":")[0]
    
        segundo = value[1].split(":")[0]
    
        if(magnitude(int(primero)) < magnitude(int(segundo))):
        
            mapping_locations_geno3[key]=value[0]
        else:
            mapping_locations_geno3[key]=value[1]
    
    return mapping_locations_geno3

## Function 8
def using_less_o_magnitude_2_dict(dictMap2):
    
    dictMap3 = {}
    
    for key in dictMap2:
        
        less_mag_Temp = dictMap2[key]
        
        mappingLocTemp = gettting_position_with_less_o_magnitude(less_mag_Temp)
        
        dictMap3[key] = mappingLocTemp
        
    return dictMap3
        
## Function 9
def making_list_correct_posi(mapping_locations_geno3):
    
    sorted_dict3 = sorted(mapping_locations_geno3.items())

    list_start_map_geno3 = []
    list_end_map_geno3 = []
    list_substraction_posi_geno = []


    for k,v in sorted_dict3:
    
        vsplit = v.split(":")
    
        start_in_geno = vsplit[0]
        end_in_geno = vsplit[1]
        
        substraction = int(end_in_geno) - int(start_in_geno)
        
        list_start_map_geno3.append(start_in_geno)
        list_end_map_geno3.append(end_in_geno)
        list_substraction_posi_geno.append(substraction)
        
    return list_start_map_geno3, list_end_map_geno3, list_substraction_posi_geno

## Function 10
def using_making_list_correct_posi_2_dict(dictMap3):
    
    dictMap4 = {}
    
    for key in dictMap3:
        
        objectFromDict = dictMap3[key]
        
        listStartTemp , listEndTemp , listSubstraction = making_list_correct_posi(objectFromDict)
        
        dictMap4[key] = (listStartTemp, listEndTemp, listSubstraction)
    
    return dictMap4

## Function 11
def getting_count_each_operon_in_dict(dict_inital):
    
    dict_operon_count = {}
    
    for key in dict_inital:
        
        operonNamesTemp = dict_inital[key][0][5]
        
        operon_count_dict = {rna:operonNamesTemp.count(rna) for rna in operonNamesTemp}
        
        operon_count_values = operon_count_dict.values()
        
        dict_operon_count[key] = operon_count_values

    return dict_operon_count

## Function 12
def dict_list_of_dfs_all_genos(dict_initial,dict_map4):
    
    # importing pandas
    import pandas as pd
    from pandas import DataFrame
    
    # dict to put dataframe
    dict_df_of_each_geno = {}
    
    # Doing loop to get all lists of each geno
    for key_initial, key_map4 in itertools.izip(dict_initial,dict_map4):
        
        #main objecto from dict initial
        mainObjectFromInitial = dict_initial[key_initial]
        
        #Creating lists of each aspect from dict_initial
        product_list = mainObjectFromInitial[0]
        start_list = mainObjectFromInitial[1]
        end_list = mainObjectFromInitial[2]
        strand_list = mainObjectFromInitial[3]
        substraction_list = mainObjectFromInitial[4]
        
        #main objecto from dict initial
        mainObjectFromMap4 = dict_map4[key_map4]
        
        #Creating lists of each aspect from dictMap4
        list_start_map_geno = mainObjectFromMap4[0]
        list_end_map_geno = mainObjectFromMap4[1]
        list_substract_geno = mainObjectFromMap4[2]
        

        #creating temp columns 
        dColumnsTemp = {

            '1_Product' : product_list,
            '2_Start': start_list,
            '3_End': end_list,
            '4_Substraction': substraction_list,
            '5_Strand_dir' : strand_list,
            '6_Start_in_geno': list_start_map_geno,
            '7_End_in_geno': list_end_map_geno,
            '8_Substract_GenoPosi': list_substract_geno

        }
        
        dfTemp = pd.DataFrame(dColumnsTemp)
        
        alert_message_lenght_operon = ""
        
        if(len(mainObjectFromInitial[5]) == len(mainObjectFromInitial[6])):

            additionalTemp = pd.DataFrame({'8_Operon_names' : mainObjectFromInitial[5],
                                       '9_NumberofSubunitsbyOperon' : mainObjectFromInitial[6]})

            newTemp = pd.concat([dfTemp,additionalTemp], axis=1)
            alert_message_lenght_operon = "It was possible to add aditional columns: OperonNames and NumberSubunits for " + key_initial 
        else:
            newTemp = dfTemp
            alert_message_lenght_operon = "It was not possible to add additional columns for " + key_initial
            
        
        key_initial_list = key_initial.split()
        
        #adding geno name
        additionalName = pd.DataFrame({'0_GenoName' : key_initial_list})
        newTemp2 = pd.concat([newTemp,additionalName],axis=1)
        
        dict_df_of_each_geno[key_initial] = newTemp2
        print alert_message_lenght_operon
        
    return dict_df_of_each_geno

## Function 13
def complete_list_dfs(dict_test_list_dfs):
    
    all_dfs = []
    for key in dict_test_list_dfs:
        
        dfTemp = dict_test_list_dfs[key]
        all_dfs.append(dfTemp)
        
    return all_dfs

## Function 14
def multiple_dfs(df_list, sheets, file_name, spaces):
    
    # It was gotten from http://stackoverflow.com/questions/32957441/
    #putting-many-python-pandas-dataframes-to-one-excel-worksheet
    
    # importing pandas
    import pandas as pd
    from pandas import DataFrame
    
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')   
    row = 0
    for dataframe in df_list:
        dataframe.to_excel(writer,sheet_name=sheets,startrow=row , startcol=0)   
        row = row + len(dataframe.index) + spaces + 1
    writer.save()

##############################################################
### Functions de Ribosomal Operons Multiple Files To Excel => Second Properties
##############################################################

## Function 15
def Processing_folder_paths(path):
    
    i = 0
    folders = []
    
    for folder in os.walk(path):

        i += 1

        if(i == 1):
            mssg = "All folder paths were obtained"
            folders = folder[1]
            break
            
    if(len(folders) == 0):
        mssg = "There are not folder in this PATH"
        
    print mssg
    return folders

## Function 16
def getting_moves_list_from_foldersPaths(folderPaths):
    
    moves_list = []
    cut_list = []
    
    for folder in folderPaths:
        
        cut_window = int(folder.split("_")[1])
        move_window = int(folder.split("_")[2])
        moves_list.append(move_window)
        cut_list.append(cut_window)
    
    return sorted(moves_list), sorted(cut_list)

## Function 17 , Auxiliar function of Function 18 => getting_positions_from_gb_file
def doing_operations_for_single_move(moves_list, startPositions, endPositions):
    
    dict_operations_for_single_move = {}
    
    for move in moves_list:
        
        moveTemp = move
        #print moveTemp
        
        listTempFrag_gbFile_start = []
        listTempFrag_gbFile_end = []
        listTemp_desv_start = []
        listTemp_desv_end = []
        
        for valueStart, valueEnd in itertools.izip(startPositions,endPositions):
            
            startTemp = int(valueStart)
            endTemp = int(valueEnd)

            fragment_in_gbFile_start = int(float(startTemp / moveTemp))
            fragment_in_gbFile_end = int(float(endTemp / moveTemp))

            # desviation to left according to start and end in gbFile
            desv_start = ( startTemp - (fragment_in_gbFile_start * moveTemp) )
            desv_end = ( endTemp - (fragment_in_gbFile_end * moveTemp) )

            listTempFrag_gbFile_start.append(fragment_in_gbFile_start)
            listTempFrag_gbFile_end.append(fragment_in_gbFile_end)
            listTemp_desv_start.append(desv_start)
            listTemp_desv_end.append(desv_end)
        
        dict_operations_for_single_move[moveTemp] = (listTempFrag_gbFile_start, listTempFrag_gbFile_end,
                                                          listTemp_desv_start, listTemp_desv_end)
        
        
    return dict_operations_for_single_move

## Function 18
def getting_positions_from_gb_file(dict_inital, moves_list):
    
    dict_positions_for_each_geno = {}
    
    for key in dict_inital:
        
        #print key
        
        startPositions = dict_inital[key][1]
        #print startPositions
        
        endPositions = dict_inital[key][2]
        #print endPositions

        dictOperationsTemp = doing_operations_for_single_move(moves_list,startPositions,endPositions)
        dict_positions_for_each_geno[key] = dictOperationsTemp
        #print "\n"
    
    mssg1 = "The keys of dict_positions_for_each_geno (Bigger Dict) are: " + ' '.join(dict_positions_for_each_geno.keys())
    mssg2 = "The keys of the Smaller Dict are: " + ''.join(str(moves_list))
    
    print mssg1
    print mssg2
    
    return dict_positions_for_each_geno

## Function 19, It can print messages to see if everything is OK
def f1_getting_definitive_operonaNames_Frag(data1,data2,data3,data4):

    finalPosiDefinitive = []

    for d1,d2,d3,d4 in itertools.izip(data1, data2, data3, data4):

        posiStart = d1
        posiEnd = d2

        desvStart = d3
        desVEnd = d4

        posiDefinitive = 0

        mssg = ""

        if(posiStart != posiEnd):
            if(desvStart < desVEnd):
                posiDefinitive = posiStart
                finalPosiDefinitive.append(posiDefinitive)
                mssg = "Everything Ok"
            else:
                posiDefinitive = posiEnd
                finalPosiDefinitive.append(posiDefinitive)
                mssg = "Everything Ok"
        elif(posiStart == posiEnd):
            posiDefinitive = posiStart
            finalPosiDefinitive.append(posiDefinitive)
            mssg = "Everything Ok"
        else:
            mssg = "There was a problem. Definitive frag could not be assigned"      

    print mssg
        #print operonNameTemp, posiDefinitive , mssg
    return finalPosiDefinitive

## Function 20

# entrada a function are the dict and the move number which are any of the following:
# [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000]
def using_f1_over_all_genos(dict_positions_for_each_geno,move_number):
    
    dict_using_f1_for_all_genos = {}
    
    for key2 in dict_positions_for_each_geno:
        
        # Getting posi info after using key2
        mainObject2 = dict_positions_for_each_geno[key2][move_number]
        
        data1 = mainObject2[0]
        data2 = mainObject2[1]
        data3 = mainObject2[2]
        data4 = mainObject2[3]
        
        
        finalPosiDefinitiveTemp = f1_getting_definitive_operonaNames_Frag(data1,data2,data3,data4)
        
        dict_using_f1_for_all_genos[key2] = finalPosiDefinitiveTemp
        
    return dict_using_f1_for_all_genos    

## Function 21
def f2_getting_operon_list_with_disposition_and_frags(data_operon_ids,startData,endData,finalPosiDefinitive):
    
    startDataNext = startData
    del(startDataNext[0])
    
    
    # Big list
    list_operons = []
    # Small list
    list_subNames = []
    
    # Big list
    list_posiDefiFrags = []
    # Small list
    list_Frags = []

    for subunit,startNext,end, posiFrag in itertools.izip_longest(data_operon_ids,startDataNext,endData,
                                                                  finalPosiDefinitive):

        subunitTemp = subunit.split(" ")[0]

        #print subunitTemp, start , end
        if(startNext != None):

            startNextTemp = int(startNext)
            endTemp = int(end)

            #print startNext , end

            substraction = startNextTemp - endTemp

            #print substraction

            list_subNames.append(subunitTemp)
            list_Frags.append(posiFrag)
            substraction = startNextTemp - endTemp

            if(substraction > 500):
                list_operons.append(list_subNames)
                list_posiDefiFrags.append(list_Frags)

                list_subNames = []
                list_Frags = []
        else:
            if(end == endData[-1]):
                list_subNames.append(subunitTemp)
                list_Frags.append(posiFrag)

                list_operons.append(list_subNames)
                list_posiDefiFrags.append(list_Frags)

                #print startNext , end
    
    finalDispositionList = []
    
    for in_list in list_operons:
        
        if(len(in_list) == 3):
            finalDisposition = in_list[0] + "-" + in_list[1] + "-" + in_list[2]
            finalDispositionList.append(finalDisposition)
        elif(len(in_list) == 4):
            finalDisposition = in_list[0] + "-" + in_list[1] + "-" + in_list[2] + "-" + in_list[3]
            finalDispositionList.append(finalDisposition)
        #else:
            #print "There was a problem giving disposition of subunits"
            #print ""

    return finalDispositionList, list_posiDefiFrags

## Function 22
def using_f2_for_all_genos(dict_initial,dict_using_f1_for_all_genos):
    
    dict_using_f2_for_all_genos = {}
    
    ## Looping over dict_initial and dict_using_f1_for_all_genos
    for key1, key2 in itertools.izip(dict_initial,dict_using_f1_for_all_genos):
        
        mainObject1 = dict_initial[key1]
        mainObject2 = dict_using_f1_for_all_genos[key2]
        
        # Getting data for function f1 from mainObject1
        data_operon_ids = mainObject1[0]

        # Getting startData and endData from mainObject1
        startData = mainObject1[1]
        endData = mainObject1[2]
        
        # Getting finalPosiDefinitive from dict_using_f1_for_all_genos
        finalPosiDefinitive = mainObject2
        
        # Using f2
        operons_f2Temp, posiDefiFrags_f2Temp = f2_getting_operon_list_with_disposition_and_frags(data_operon_ids,
                                                                                                  startData,
                                                                                                  endData,
                                                                                                  finalPosiDefinitive)
        
        # Getting data for f3 function
        count_operon_namesTemp = mainObject1[5]
        count_operon_listTemp = mainObject1[6]
        
        # GenoName
        genoNameTemp = str(key1)
        
        if( len(operons_f2Temp) == 0 or len(operons_f2Temp) != len(posiDefiFrags_f2Temp) ):
            print "There is a problem with dispositions in " + key1
            
            # Inserting lists to dict
            dict_using_f2_for_all_genos[key1] = (operons_f2Temp,posiDefiFrags_f2Temp,count_operon_namesTemp,
                                                    count_operon_listTemp,genoNameTemp)
        else:
            
            # Inserting lists to dict
            dict_using_f2_for_all_genos[key1] = (operons_f2Temp,posiDefiFrags_f2Temp,count_operon_namesTemp,
                                                    count_operon_listTemp,genoNameTemp)
    
    return dict_using_f2_for_all_genos

## Function 22 Version 2
def using_f2_for_all_genos_v2(dict_initial,dict_using_f1_for_all_genos):
    
    dict_using_f2_for_all_genos = {}

    set_dict_initial = set(dict_initial)
    set_dict_f1 = set(dict_using_f1_for_all_genos)
    
    ## Looping over intersection dict_initial and dict_using_f1_for_all_genos
    for key in set_dict_initial.intersection(set_dict_f1):
        
            mainObject1 = dict_initial[key]
            mainObject2 = dict_using_f1_for_all_genos[key]

            # Getting data for function f1 from mainObject1
            data_operon_ids = mainObject1[0]

            # Getting startData and endData from mainObject1
            startData = mainObject1[1]
            endData = mainObject1[2]

            # Getting finalPosiDefinitive from dict_using_f1_for_all_genos
            finalPosiDefinitive = mainObject2

            # Using f2
            operons_f2Temp, posiDefiFrags_f2Temp = f2_getting_operon_list_with_disposition_and_frags(data_operon_ids,
                                                                                                      startData,
                                                                                                      endData,
                                                                                                      finalPosiDefinitive)

            # Getting data for f3 function
            count_operon_namesTemp = mainObject1[5]
            count_operon_listTemp = mainObject1[6]

            # GenoName
            genoNameTemp = str(key)

            if( len(operons_f2Temp) == 0 or len(operons_f2Temp) != len(posiDefiFrags_f2Temp) ):
                print "There is a problem with dispositions in " + key

                # Inserting lists to dict
                dict_using_f2_for_all_genos[key] = (operons_f2Temp,posiDefiFrags_f2Temp,count_operon_namesTemp,
                                                        count_operon_listTemp,genoNameTemp)
            else:

                # Inserting lists to dict
                dict_using_f2_for_all_genos[key] = (operons_f2Temp,posiDefiFrags_f2Temp,count_operon_namesTemp,
                                                        count_operon_listTemp,genoNameTemp)
            

    return dict_using_f2_for_all_genos

## Function 23
def f3_final_posiFrag_mean_for_each_geno(list_posiDefiFrags_f2):
    
    ## getting numpy
    import numpy as np
    
    ## list of definitive frag position after mean operation
    final_posiFrag_mean_for_each_geno = []
    operonNumberList = []
    
    count = 0
    for operonSubUnit in list_posiDefiFrags_f2:
        
        posiFinalAfterMean = int(np.mean(operonSubUnit))
        
        final_posiFrag_mean_for_each_geno.append(posiFinalAfterMean)
        
        count += 1
        operonNumberList.append(count)
        
    return operonNumberList, final_posiFrag_mean_for_each_geno

## Function 24
def using_f3_for_all_genos(dict_using_f2_for_all_genos):
    
    dict_using_f3_for_all_genos = {}
    
    ## Looping over dict_using_f2_for_all_genos
    for key in dict_using_f2_for_all_genos:
        
        mainObject = dict_using_f2_for_all_genos[key]
        
        list_posiDefiFrags_f2Temp = mainObject[1]
        
        # using f3 over argument that is input
        operNumberListTemp, finalPosiFrag_meanTemp = f3_final_posiFrag_mean_for_each_geno(list_posiDefiFrags_f2Temp)
        
        dict_using_f3_for_all_genos[key] = (operNumberListTemp,finalPosiFrag_meanTemp)
        
    return dict_using_f3_for_all_genos

## Function 25
def f4_doing_data_frame_for_each_geno(operonNumberList,list_operons_f2,
                                      final_posiFrag_f2,count_operon_names,
                                      count_operon_list,genoName):
    
     # importing pandas
    import pandas as pd
    from pandas import DataFrame
    
    len1 = len(operonNumberList)
    len2 = len(list_operons_f2)
    len3 = len(final_posiFrag_f2)
    
    if(len1 != len2 and len1 == len3):
        
        dColumns={

            '1_OperoNumbers' : operonNumberList,
            '2_Posi_Definitive_in_Image': final_posiFrag_f2

        }
        
        df = pd.DataFrame(dColumns)
        
        additionalDisposition = pd.DataFrame({'3_DispositionOperons' : list_operons_f2})
        newTemp1 = pd.concat([df,additionalDisposition], axis=1)
        
        # Adding other things
        additionalTemp = pd.DataFrame({'4_Operon_names' : count_operon_names,
                                           '5_NumberofSubunitsbyOperon' : count_operon_list})
        newTemp2 = pd.concat([newTemp1,additionalTemp], axis=1)

        #adding geno name
        additionalName = pd.DataFrame({'6_GenoName' : genoName},index=[0])
        newTemp3 = pd.concat([newTemp2,additionalName], axis=1) 
        
        dataFrameToReturn = newTemp3

    else:
        dColumns={

            '1_OperoNumbers' : operonNumberList,
            '2_Posi_Definitive_in_Image': final_posiFrag_f2,
            '3_DispositionOperons': list_operons_f2
        }

        df = pd.DataFrame(dColumns)
        
        # Adding other things
        additionalTemp = pd.DataFrame({'4_Operon_names' : count_operon_names,
                                           '5_NumberofSubunitsbyOperon' : count_operon_list})
        newTemp1 = pd.concat([df,additionalTemp], axis=1)

        #adding geno name
        additionalName = pd.DataFrame({'6_GenoID' : genoName},index=[0])
        newTemp2 = pd.concat([newTemp1,additionalName], axis=1)
        
        dataFrameToReturn = newTemp2

    return dataFrameToReturn

## Function 26
def using_f4_for_all_genos(dict_using_f2_for_all_genos,dict_using_f3_for_all_genos):
    
    ## Generating dicts of dataframes
    dict_using_f4_for_all_genos = {}
    
    ## Looping over dict_using_f2_for_all_genos and dict_using_f3_for_all_genos
    for key1,key2 in itertools.izip(dict_using_f2_for_all_genos,dict_using_f3_for_all_genos):
        
        mainObject1 = dict_using_f2_for_all_genos[key1]
        mainObject2 = dict_using_f3_for_all_genos[key2]
        
        # Getting data that f4 received as input
        operDisposition_f2 = mainObject1[0]
        count_operon_names = mainObject1[2]
        count_operon_list = mainObject1[3]
        genoName = mainObject1[4]
        
        operonNumberList_f3 = mainObject2[0]
        operFinalPosiAfterMean_f3 = mainObject2[1] 
        
        dfTemp = f4_doing_data_frame_for_each_geno(operonNumberList_f3,operDisposition_f2,
                                                  operFinalPosiAfterMean_f3,count_operon_names,
                                                  count_operon_list,genoName)
        
        dict_using_f4_for_all_genos[key1] = dfTemp
    
    return dict_using_f4_for_all_genos



    