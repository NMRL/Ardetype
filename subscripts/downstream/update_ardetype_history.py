import sys, pandas as pd
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to scripts
report_time = uu.get_current_timestamp()
arg_dict = {
    "ardetype"    : ["--ard", 'Full path to ardetype report (e.g. ardetype_report.csv)']
}

proc_dict = uu.proc_dict
proc_dict['ard'] = uu.proc_dict['default'].copy()

##########
#Functions
##########

def custom_ardetype(df:pd.DataFrame):
    '''Custom function to perform method-based result aggregation from ardetype report'''
    df = df.applymap(str)
    tools = [col.replace("method|", "") for col in df.columns if "method" in col]
    if tools: #if species-specific typing was performed
        df["methods"] = df[[col for col in df.columns if "method|" in col]].agg(' ; '.join, axis=1) #CONCATENATE METHODS
        df["types"] = df[[col for col in df.columns if "type|" in col]].agg(' ; '.join, axis=1) #CONCATENATE TYPES
        df["types"] = [(' ; ').join([tools[i]+": "+value for i,value in enumerate(row.split(' ; '))]) for row in df['types']] #ADD TOOL-BASED ANNOTATION TO EACH TYPE
        df["references"] = df[[col for col in df.columns if "reference|" in col]].agg(' ; '.join, axis=1) #CONCATENATE REFERENCES
        df.drop(list(df.filter(regex = '(method\||type\||reference\|)')), axis = 1, inplace = True) #REMOVE UNUSED COLUMNS
        df = df.replace('((; |)nan( |)|; [A-z0-9]*: nan |[A-z0-9]*: nan ; |[\-A-z0-9]*: nan)','', regex=True) #REMOVE REDUNDANT INFORMATION
        df = df.replace('(; $|^; )','', regex=True) #remove redundant separators
        df.drop('taxid', inplace=True, axis=1) #drop unused column
    else:
        #blank entries for samples where no specific data is available
        df['methods'] = ["" for _ in df.index]
        df['types'] = ["" for _ in df.index] 
        df['references'] = ["" for _ in df.index]
        df.fillna("", inplace=True)
    return df

##############
#Runtime logic
##############

if __name__ == '__main__':
    parser = uu.parse_arguments(arg_dict)
    proc_dict['ard']['custom_functions'] = [custom_ardetype]
    uu.create_backup(full_path)
    current_tables = uu.find_current_tables(full_path, arg_dict)
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
