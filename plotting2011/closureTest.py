
import sys,os, math
import subprocess as sp

def getDict ( file_name, data_name ) :
    file = open ( file_name, "r" ) 
    file_info = file.read() 
    data_info = file_info.split(data_name)[1].split("\n\n")[0]
    
    list = []
    dict = {} 

    data_info_lines = data_info.split("\n")
    for line in data_info_lines : 
        split_line = line.split()
        if len (split_line) == 0: continue
        if split_line[0] == "variableName" : continue
        
        variable = split_line [0]
        npass    = split_line [5]
        eNpass   = split_line [6] 
        list.append ( variable ) 
        dict [ variable ] = [ npass, eNpass ]
    return list, dict 
        

cej_file_name = "/afs/cern.ch/user/e/eberry/scratch0/LQDATA//eejj_analysis/QCDClosureTest_cej_DataAndMC/output_cutTable_lq_QCDClosureTest_cej_DataAndMC/analysisClass_lq_eejj_QCDClosureTest_tables.dat"

ccj_file_name = "/afs/cern.ch/user/e/eberry/scratch0/LQDATA//eejj_analysis/QCDClosureTest_ccj_DataOnly/output_cutTable_lq_QCDClosureTest_ccj_DataOnly/analysisClass_lq_eejj_QCDClosureTest_tables.dat"

ccjp_file_name = "/afs/cern.ch/user/e/eberry/scratch0/LQDATA//eejj_analysis/QCDClosureTest_ccj_DataOnly_Plus1Sigma/output_cutTable_lq_QCDClosureTest_ccj_DataOnly/analysisClass_lq_eejj_QCDClosureTest_tables.dat"

ccjm_file_name = "/afs/cern.ch/user/e/eberry/scratch0/LQDATA//eejj_analysis/QCDClosureTest_ccj_DataOnly_Minus1Sigma/output_cutTable_lq_QCDClosureTest_ccj_DataOnly/analysisClass_lq_eejj_QCDClosureTest_tables.dat"


variable_list , cej_allbkg_dict = getDict ( cej_file_name , "ALLBKG" ) 
variable_list , cej_data_dict   = getDict ( cej_file_name , "DATA"   ) 
variable_list , ccj_data_dict   = getDict ( ccj_file_name , "DATA"   ) 
variable_list , ccjp_data_dict  = getDict ( ccjp_file_name, "DATA"   ) 
variable_list , ccjm_data_dict  = getDict ( ccjm_file_name, "DATA"   ) 

tab_length = 7

for variable in variable_list[6:] : 
    
    predicted        = ccj_data_dict   [variable][0]

    predicted_p      = ccjp_data_dict  [variable][0]
    predicted_m      = ccjm_data_dict  [variable][0]

    ePredicted = max ( abs ( float (predicted) - float(predicted_m) ) , 
                       abs ( float (predicted) - float(predicted_p) ) ) 
                                       
    observed_data    = cej_data_dict   [variable][0]
    eObserved_data   = cej_data_dict   [variable][1]
                                       
    observed_allbkg  = cej_allbkg_dict [variable][0]
    eObserved_allbkg = cej_allbkg_dict [variable][1]

    observed         = float(observed_data) - (float(observed_allbkg) / 1000.) 
    eObserved        = math.sqrt ( ( (float(eObserved_data  )      ) * ( float(eObserved_data  )        ) ) + 
                                   ( (float(eObserved_allbkg)/1000.) * ( float(eObserved_allbkg) / 1000.) ) )
    
    ratio  = float ( predicted ) / float ( observed ) 
    eRatio = ratio * math.sqrt ( ( ( float(ePredicted) / float(predicted ) ) * ( float(ePredicted) / float(predicted ) ) ) + 
                                 ( ( float(eObserved ) / float(observed  ) ) * ( float(eObserved ) / float(observed  ) ) ) )

    space = "\t\t\t"
    if len (variable.strip() ) < tab_length:     space = space + "\t"
    if len (variable.strip() ) > tab_length * 2: space = space[:-1]
    if len (variable.strip() ) > tab_length * 3: space = space[:-1]

    print variable.strip(), space , predicted, "\t+/-\t",  ePredicted, "\t", observed, "\t +/- \t",  eObserved , "\t", ratio, "\t +/- \t",  eRatio


