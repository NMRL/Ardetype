# Short illustration of use cases for update_utilities.py

## Scope
***update_utilities.py*** module contains functions and variable templates that can be useful to aggregate tabular files (tsv, csv etc.) into single file to be loaded into relational database.

## Basic use case (with default configuration)
- Simply apply the functions from update_utilities.py
```
import sys, pandas as pd
sys.path.insert(0,'/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/')
from subscripts.downstream import update_utilities as uu

#################
#Global variables
#################

full_path = uu.get_folder_path(__file__) #path to the specific script
report_time = uu.get_current_timestamp() #string timestamp

#defining command-line arguments
arg_dict = {
    "<report_prefix>"    : ["--<abbreviation>", "Example help message specifying the expected file name, e.g. <report_prefix>_report.csv"],
}

##############
#Runtime logic
##############

if __name__ == '__main__':
    #parse command-line arguments
    parser = uu.parse_arguments(arg_dict)
    #create backup folder in the folder where the script is located (if does not exist)
    uu.create_backup(full_path)
    #detect current aggregated files (must exist) in the folder where the script is located
    current_tables = uu.find_current_tables(full_path, arg_dict)
    #perform backups and add information from new files to aggregated files
    uu.update_files(arg_dict, parser, current_tables, report_time, full_path, proc_dict)
```

## Adjusting default configuration to fit specific case
```
...
#get template proc_dict from update_utilities module
proc_dict = uu.proc_dict

#add non-default configuration as copy of default 
proc_dict['hamr'] = proc_dict['default'].copy() #copy required to avoid overwriting the default

#changing the default delimiter
proc_dict['hamr']["delimiter"] = '\t' #to read table correctly

#change default sample identifier column
proc_dict['hamr']["id_column"] = 'input_file_name'

#add more substrings to remove in id_column
proc_dict['hamr']["id_regex"].append('.amr.alignment')

#change default primary key
proc_dict['hamr']["primary_key"] = None #if set to None, all columns will be used to infer uniqueness of rows

#Change default history file extenstion
proc_dict['hamr']["history_ext"] = '.tsv'

#loop can be used if several inputs require similar processing settings
for report in ['rf', 'pf', 'afm']:
    proc_dict[report] = proc_dict['default'].copy()
    proc_dict[report]['primary_key'] = None
...
```
- See ***update_resistance_history.py*** for full example

## Applying custom functions
```
...
#get template proc_dict from update_utilities module
proc_dict = uu.proc_dict
#add non-default configuration as copy of default 
proc_dict['ard'] = uu.proc_dict['default'].copy()

##########
#Functions
##########

#defining custom function to extend default dataframe processing
def custom_ardetype(df:pd.DataFrame):
    '''Custom function that takes in dataframe and returns it after processing'''
    #custom dataframe processing logic (hopefully preserving the primary keys)
    return df

##############
#Runtime logic
##############

if __name__ == '__main__':
    #parse command-line arguments
    parser = uu.parse_arguments(arg_dict)
    #Adding list of custom function for specific input file
    proc_dict['ard']['custom_functions'] = [custom_ardetype]
    ...
```
- See ***update_ardetype_history.py*** for full example
