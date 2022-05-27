import sys, argparse, re, logging, urllib, pandas as pd,  plotly.graph_objects as go, json, os
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from igraph import Graph, EdgeSeq


class Interface:
    '''Class is set to represent command line interface to interact with local database'''
    def __init__(self, alphabet, valid_input_type, accession_pattern):
        self.alphabet = alphabet
        self.valid_input_type = valid_input_type
        self.accession_pattern = accession_pattern
        
    def parse_arguments(self):
        '''Method is used to parse command-line arguments, allowing to interact with the database from script.'''
        parser = argparse.ArgumentParser(description='This script can be used to work with local database of bacterial reference genomes in fasta format.')
        req_arg_grp = parser.add_argument_group('Required arguments')
        arg_dict = {
            "--add-fasta": f"Add new record(s) from locally stored fasta file & metadata file of {self.valid_input_type} format. Expected: file_name.{self.valid_input_type}",
            "--add-ncbi": f"Add new record from NCBI RefSeq database by correct accession number. Expected: valid_accession_number",
            "--add-ncbi-list": f"Add new record(s) from NCBI RefSeq database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format (single column with one header row). Expected: path_to_{self.valid_input_type}",
            "--exp-fasta": f"Create a fasta file containing sequences from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--exp-meta": f"Create a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--exp-records": f"Create a fasta file containing sequences & a csv file containing taxonomic information from local database based on a list of correct accession numbers provided as a file of {self.valid_input_type} format. Expected: path_to_{self.valid_input_type}",
            "--rm-record": f"Remove records from local database if given correct accession number. Expected: valid_accession_number_from_local_db",
            "--ch-header": f"Replace an accession number that exists in local database with user-provided accession number if provided accession number was found in local database. Expected: valid_accession_from_local_db,new_accession",
            "--ch-tax": f"Replace a taxonomy string that exists in local database with user-provided taxonomy string if provided accession number was found in local database. Expected: valid_accession_from_local_db,new_taxonomy_string",
            "--view-data": "Visualize the contents of the local database as a tree chart, showing the number of records that belong to each taxonomic group."
        }

        for arg in arg_dict.keys():
            if arg != "--view-data":
                req_arg_grp.add_argument(arg, metavar='\b', help = arg_dict[arg],required=False)  
            else:
                req_arg_grp.add_argument(arg, help = arg_dict[arg], action='store_true', required=False) #To allow flag instead of cmd argument

        if len(sys.argv)==1: #If no arguments were provided - print help & exit
            parser.print_help(sys.stderr)
            sys.exit(1)
        self.args = parser.parse_args()
        return self.args

    def check_format(self, raw_accession):
        '''Method checks if provided accession number matches RefSeq accession number conventions.'''
        if re.match(self.accession_pattern, raw_accession): return True
        else: return False

    def check_alphabet(self, raw_sequence):
        '''Method checks if sequence that corresponds to the id contains only allowed characters.'''
        for letter in str(raw_sequence):
            if letter not in self.alphabet:
                return False
        return True

    def check_file_type(self, file_name):
        '''Method checks if provided file type matches allowed file type.'''
        return all([os.path.exists(f'{file_name}.{self.valid_input_type}'), os.path.exists(f'{file_name}.fasta')])

    def read_local(self, file_name):
        '''Method is used to read local sequence and taxonomy information.'''
        local_seq = SeqIO.to_dict(SeqIO.parse(f'{file_name}.fasta', "fasta"))
        local_tax_df = pd.read_csv(f'{file_name}.{self.valid_input_type}', header=[0])
        return local_seq, local_tax_df


class Database:
    '''Class is set to represent local database of sequences in fasta format along with annotation csv table'''
    def __init__(self, db_name):
        '''Constructor for the database class'''
        self.sequence_file = Path(f'{db_name}.fasta')
        self.taxonomy_file = Path(f'{db_name}.csv')
        
    def create_db_files(self):
        '''Method is used to create database file in the directory where the script is.'''
        try:
            self.sequence_file.touch(exist_ok=False) #Create db file if it does not yet exist, else raise an exception
            self.taxonomy_file.touch(exist_ok=False)
            tax_file_columns = pd.DataFrame(columns=["accession_number","taxonomy_string"]) #Initialize empty dataframe
            tax_file_columns.to_csv(self.taxonomy_file, index=False, header=True) #Store initialized dataframe in a file (keep headers)
        except:
            print(f'Current database files: {self.sequence_file}; {self.taxonomy_file}') #This message is printed if database files already exist

    def add_record(self, header, sequence, taxonomy, description):
        '''Method is used to add new record to database files (taxonomy file & fasta file).'''
        new_record_seq = SeqIO.SeqRecord(Seq(sequence),id=header, description=description) #Initialize new SeqRecord class object to store sequence information untill it is added to the file 
        new_record_tax = pd.DataFrame.from_dict({"accession_number":[header],"taxonomy_string":[taxonomy]}) #Initialize a dataframe and add accession number and taxonomic of the added sequence to it
        with open(self.sequence_file, "a+") as seq_db: #Open database file for sequences
            SeqIO.write(new_record_seq, seq_db, "fasta") #Add new sequence record to the file & update the sequence file content
        tax_db = pd.read_csv(self.taxonomy_file, header=[0]) #Open database file for taxonomy
        tax_db = tax_db.append(new_record_tax) #Add new taxonomy record to the file
        tax_db.to_csv(self.taxonomy_file, index=False) #Update the taxonomy file content
        
    def calculate_content(self, taxonomy_string=None, reduce=False, new_accession=None, accession=None, recalc_accessions=False):
        '''Method is used to calculate current database content by taxonomic group and store in a file.'''
        if recalc_accessions:
            with open(f'accession_map.json', 'r+') as f4:
                acc_dict = json.load(f4)
            with open(f'encode_dict.json', "r+") as f3:
                encode_dict = json.load(f3)
            taxonomy_list = taxonomy_string.split("|")
            for taxon in taxonomy_list:
                if accession in acc_dict[str(encode_dict[taxon])]:
                    acc_dict[str(encode_dict[taxon])].remove(accession)
                    acc_dict[str(encode_dict[taxon])].append(new_accession)
            with open(f'accession_map.json', "w+") as file:
                json.dump(acc_dict, file)
            return acc_dict
        try:
            with open(f'adj_set.json', "r+") as f1:
                adj_set = json.load(f1)
            with open(f'count_dict.json', "r+") as f2:
                count_dict = json.load(f2)
            with open(f'encode_dict.json', "r+") as f3:
                encode_dict = json.load(f3)
            with open(f'accession_map.json', 'r+') as f4:
                acc_dict = json.load(f4)
            adj_set = set(tuple(pair) for pair in adj_set)
            if taxonomy_string:
                if not reduce:
                    taxonomy_list = taxonomy_string.split("|")
                    new_adj_dict = {taxonomy_list[i]:taxonomy_list[i+1] for i in range(len(taxonomy_list)-1)}
                    for taxon in taxonomy_list:
                        try:
                            count_dict[taxon] += 1
                            encode_dict[taxon]
                        except KeyError:
                            encode_dict[taxon] = max(encode_dict.values()) + 1
                            count_dict[taxon] = 1
                    for parent,child in new_adj_dict.items():
                        adj_set.add((encode_dict[parent],encode_dict[child]))
                    
                    tax_db = pd.read_csv(self.taxonomy_file, header=[0])
                    new_split_product = tax_db['taxonomy_string'].str.split("|",expand=True)
                    new_split_product['accession_number'] = tax_db['accession_number']

                    for taxon in taxonomy_list:
                        new_accession = list(new_split_product[new_split_product.isin([taxon]).any(axis=1)]['accession_number'])
                        if taxon not in acc_dict.keys():
                            acc_dict[encode_dict[taxon]] = new_accession
                        else:
                            acc_dict[encode_dict[taxon]] = list(set(acc_dict[taxon] + new_accession))

                else:
                    os.remove('adj_set.json')
                    os.remove('count_dict.json')
                    os.remove('encode_dict.json')
                    os.remove('accession_map.json')
                    raise FileNotFoundError
                        
                with open(f'adj_set.json', "w+") as file:
                        json.dump(list(adj_set), file)
                with open(f'count_dict.json', "w+") as file:
                    json.dump(count_dict, file)
                with open(f'encode_dict.json', "w+") as file:
                    json.dump(encode_dict, file)
                with open(f'accession_map.json', "w+") as file:
                    json.dump(acc_dict, file)
            return adj_set, count_dict, acc_dict, encode_dict
        except FileNotFoundError:
            tax_db = pd.read_csv(self.taxonomy_file, header=[0])
            split_product = tax_db['taxonomy_string'].str.split("|",expand=True)
            new_split_product = tax_db['taxonomy_string'].str.split("|",expand=True)
            new_split_product['accession_number'] = tax_db['accession_number']
            adj_set = set()
            count_dict = {}
            for col in split_product.columns:
                count_dict.update(dict(split_product[col].value_counts()))
            
            encode_dict = {node:i for i,node in enumerate(count_dict.keys())}
            
            acc_dict = {}
            for taxon in encode_dict.keys():
                new_accession = list(new_split_product[new_split_product.isin([taxon]).any(axis=1)]['accession_number'])
                if taxon not in acc_dict.keys():
                    acc_dict[taxon] = new_accession
                else:
                    acc_dict[taxon] = list(set(acc_dict[taxon] + new_accession))
            
            new_acc_dict = {}
            for taxon in acc_dict.keys():
                new_acc_dict[encode_dict[taxon]] = acc_dict[taxon]

            count_dict = {key:int(value) for key,value in count_dict.items()}
            split_product = split_product.drop_duplicates()
            for i in range(len(split_product.columns)):
                for j in range(len(split_product.iloc[:,i])):
                    if i+1 != len(split_product.columns):
                        if split_product.iloc[j,i] is None or split_product.iloc[j,i+1] is None:
                            continue
                        else:
                            adj_set.add((encode_dict[split_product.iloc[j,i]], encode_dict[split_product.iloc[j,i+1]]))
            with open(f'adj_set.json', "w+") as file:
                json.dump(list(adj_set), file)
            with open(f'count_dict.json', "w+") as file:
                json.dump(count_dict, file)
            with open(f'encode_dict.json', "w+") as file:
                json.dump(encode_dict, file)
            with open(f'accession_map.json', 'w+') as file:
                json.dump(new_acc_dict, file)
            return adj_set, count_dict, acc_dict, encode_dict
        
    def find_id(self, valid_accession):
        '''Method is used to check if given id exists in the database.'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        if valid_accession in list(tax_db["accession_number"]): return True
        else: return False

    def find_tax(self, valid_accession):
        '''Method is used to get taxonomy information of a record in the database.'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]
        return tax_db.iloc[accession_ind].loc['taxonomy_string']

    def rm_record(self, valid_accession):
        '''Method is used to remove a record from local database based on accession number.'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        tax_db = tax_db[tax_db["accession_number"] != valid_accession]
        del seq_db[valid_accession]
        with open(self.sequence_file, "w+") as seq_db_file:
            SeqIO.write(seq_db.values(), seq_db_file, "fasta")
        tax_db.to_csv(self.taxonomy_file, index=False)

    def write_tax(self, valid_accession, new_taxonomy):
        '''Method is used to replace a taxonomy information of a record in the database.'''
        if "|" in new_taxonomy:
            tax_db = pd.read_csv(self.taxonomy_file, header=[0])
            accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]

            tax_db.at[accession_ind,'taxonomy_string']=new_taxonomy

        else:
            sys.exit('Invalid taxonomy string. Expected: A|B|...|C with A,B,C as taxons (length irrelevant)')
        
        tax_db.to_csv(self.taxonomy_file, index=False)

    def write_id(self, old_accession, new_accession):
        '''Method is used to replace an accession of a record in the database.'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        accession_ind = tax_db.index[tax_db["accession_number"] == old_accession].to_list()[0]
        tax_db.at[accession_ind,'accession_number']=new_accession
        seq_db[new_accession] = seq_db[old_accession]
        seq_db[new_accession].id = new_accession
        seq_db[new_accession].description = f'{new_accession} replaced manually'
        del seq_db[old_accession]
        seq_db[new_accession].id = new_accession
        with open(self.sequence_file, "w+") as seq_db_file:
            SeqIO.write(seq_db.values(), seq_db_file, "fasta")
        tax_db.to_csv(self.taxonomy_file, index=False)

    def export_fasta(self, valid_accession):
        '''Method is used to export single sequence in a form of SeqIO record.'''
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        requested_seq = seq_db[valid_accession]
        return requested_seq
               
    def export_tax(self, valid_accession):
        '''Method is used to export single taxonomy record in a form of numpy.Series'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]
        requested_tax = tax_db.iloc[accession_ind]
        return requested_tax

    def export_record(self, valid_accession):
        '''Method is used to export single sequence in a form of SeqIO record and corresponding taxonomy record in a form of numpy.Series'''
        tax_db = pd.read_csv(self.taxonomy_file, header=[0])
        seq_db = SeqIO.to_dict(SeqIO.parse(self.sequence_file, "fasta"))
        accession_ind = tax_db.index[tax_db["accession_number"] == valid_accession].to_list()[0]
        requested_tax = tax_db.iloc[accession_ind]
        requested_seq = seq_db[valid_accession]
        return requested_tax, requested_seq


class Query_ncbi:
    '''Class is set to contain methods to send requests to ncbi nucleotide database.'''
    def __init__(self, database, email, alphabet):
        self.database = database
        self.email = email
        self.alphabet = alphabet

    def check_output(self, valid_accession):
        '''Method checks if record exists in NCBI database given accession number.'''
        Entrez.email = self.email
        try:
            handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="acc") #Attempts to open connection to NCBI database
            data = handle.read() #Reads up-to-date accession from NCBI
            handle.close() #Closes connection
            if valid_accession in data: #Checks if found accession contains query
                return True
            else: 
                return False
        except urllib.error.HTTPError as e: #If query is invalid, the exception is thrown
            if str(e) == "HTTP Error 400: Bad Request":
                return False

    def get_fasta(self, valid_accession):
        '''Method attempts to get fasta sequence (header + sequence separately) given accession number.'''
        #generate request
        Entrez.email = self.email
        handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="fasta") #Attempts to open connection to NCBI database
        record = SeqIO.read(handle, "fasta") #Reads sequence in fasta format from ncbi
        handle.close() #Closes connection
        return record.seq, record.id

    def get_taxonomy(self, valid_accession):
        '''Method attempts to get taxonomy information(list) from NCBI database given accession number.'''
        Entrez.email = self.email
        handle = Entrez.efetch(db=self.database, id = valid_accession, rettype="gb") #Attempts to open connection to NCBI database
        tax_data = SeqIO.read(handle, format='genbank').annotations['taxonomy']
        handle.close() #Closes connection
        return tax_data


class Logger:
    '''Class is set to provide logging options to store history of edits done to the local database via command line'''
    def __init__(self, operation_type, valid_accession=None,  old_accession=None,  old_taxonomy=None,  new_accession=None,  new_taxonomy=None, log_file='PD2.log'):
        self.operation_type = operation_type
        self.log_file = log_file
        self._valid_accession = valid_accession
        self._old_accession = old_accession
        self._old_taxonomy = old_taxonomy
        self._new_accession = new_accession
        self._new_taxonomy = new_taxonomy
        self.log_dict = {
        "add_fasta":f"New record was added from fasta file (id: {self._valid_accession})",
        "add_ncbi":f"New record was added from ncbi database (id: {self._valid_accession})",
        "rm_record":f"Record was removed (id: {self._valid_accession})",
        "ch_header":f"ID of the record was changed from {self._old_accession} to {self._new_accession}",
        "ch_tax":f"Taxonomy information for {self._valid_accession} was changed from {self._old_taxonomy} to {self._new_taxonomy}"
    }

    def update_log_dict(self):
        self.log_dict = {"add_fasta":f"New record was added from fasta file (id: {self._valid_accession})",
        "add_ncbi":f"New record was added from ncbi database (id: {self._valid_accession})",
        "rm_record":f"Record was removed (id: {self._valid_accession})",
        "ch_header":f"ID of the record was changed from {self._old_accession} to {self._new_accession}",
        "ch_tax":f"Taxonomy information for {self._valid_accession} was changed from {self._old_taxonomy} to {self._new_taxonomy}"}

    def create_logger(self):
        '''Method is used to create and configure logger''' 
        #Create log file if it does not exist
        log_file = Path(self.log_file)
        log_file.touch(exist_ok=True)
        
        #Log changes - full documentation link - https://docs.python.org/3/howto/logging.html
        self.pd2_logger = logging.getLogger(self.operation_type)
        self.pd2_logger.setLevel(logging.DEBUG)
        self.log_file_handler = logging.FileHandler(filename=self.log_file)
        self.log_file_formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S')
        self.log_file_handler.setFormatter(self.log_file_formatter)
        self.pd2_logger.addHandler(self.log_file_handler)

    def log_change(self):
        '''Method is used to log local changes to database files'''
        self.pd2_logger.info(self.log_dict[self.operation_type])

    def close_logger(self):
        '''Method is used to close file handler after logging is done'''
        self.log_file_handler.close() #Close file after logging is complete


class Plotter:
    '''Class is set to provide visual view options for the content of the local database'''
    def __init__(self, adj_set, count_dict, accession_map, encode_dict):
        self.adj_set = adj_set
        self.count_dict = count_dict
        self.accession_map = accession_map
        self.encode_dict = encode_dict

    def display(self):
        '''Method is used to display current content of the local database'''
        #Adapted from https://plotly.com/python/tree-plots/
        nr_vertices = len(self.count_dict.keys())
        nl = '<br>'
        v_label = [f'{key}: {self.count_dict[key]}{nl}{nl.join(self.accession_map[str(self.encode_dict[key])])}' for key in self.count_dict.keys()]
        ball_names = [f'{key}: {self.count_dict[key]}' for key in self.count_dict.keys()]
        graph = Graph()
        graph.add_vertices(nr_vertices)
        graph.add_edges(list(self.adj_set))
        lay = graph.layout('tree')
        max_y = max([lay[k][1] for k in range(nr_vertices)])

        es = EdgeSeq(graph) # sequence of edges
        E = [e.tuple for e in graph.es] # list of edges

        Xn = [lay[k][0] for k in range(nr_vertices)]
        Yn = [2*max_y-lay[k][1] for k in range(nr_vertices)]
        Xe = []
        Ye = []
        for edge in E:
            Xe+=[lay[edge[0]][0],lay[edge[1]][0], None]
            Ye+=[2*max_y-lay[edge[0]][1],2*max_y-lay[edge[1]][1], None]

        labels = v_label
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=Xe, y=Ye, mode='lines', line=dict(color='rgb(210,210,210)', width=1), hoverinfo='none'))
        fig.add_trace(go.Scatter(x=Xn, y=Yn, mode='markers', name='bla', marker=dict(symbol='circle-dot', size=100, color='#6175c1',    #'#DB4551',
                                        line=dict(color='rgb(50,50,50)', width=1)), text=v_label, hoverinfo='text', opacity=0.8))


        def make_annotations(pos, text, font_size=15, font_color='rgb(0,0,0)'):
            L=len(pos)
            if len(text)!=L:
                raise ValueError('The lists pos and text must have the same len')
            annotations = []
            for k in range(L):
                annotations.append(
                    # or replace labels with a different list for the text within the circle
                    dict(text=ball_names[k], x=pos[k][0], y=2*max_y-lay[k][1], xref='x1', yref='y1', font=dict(color=font_color, size=font_size), showarrow=False)
                )
            return annotations

        axis = dict(showline=False, zeroline=False, showgrid=False, showticklabels=False,) # hide axis line, grid, ticklabels and  title

        fig.update_layout(title= 'Local database summary', annotations=make_annotations(lay, v_label), font_size=14, showlegend=False,
                      xaxis=axis, yaxis=axis, margin=dict(l=40, r=40, b=85, t=100), hovermode='closest', plot_bgcolor='rgb(143, 150, 196)')
        fig.show()


if __name__ == "__main__":
    my_interface = Interface('ACGT', 'csv', r"^[A-z]{2}_[A-Z0-9]*$")
    arg_dict = vars(my_interface.parse_arguments())
    my_logger = Logger(None, valid_accession="")
    my_logger.create_logger()
    if any([True if arg else False for arg in arg_dict.values()]):
        my_database = Database('my_database')
        my_database.create_db_files()
        my_qry = Query_ncbi('nucleotide', 'email@example.com', 'ACGT')

    if arg_dict['add_fasta']:
        file_check = my_interface.check_file_type(f'{arg_dict["add_fasta"]}')
        if file_check:
            records = my_interface.read_local(arg_dict['add_fasta'])
            sequences, taxonomy = records[0], records[1]
            for key in sequences.keys():
                acc_nr = sequences[key].id
                sqnc = sequences[key].seq
                tax_ind = taxonomy.index[taxonomy['accession_number'] == key].to_list()[0]
                tax_str = taxonomy.iloc[tax_ind]['taxonomy_string']
                ab_check = my_interface.check_alphabet(sqnc)
                id_check = my_interface.check_format(acc_nr)
                dupl_check = my_database.find_id(acc_nr)
                descr = f'from local file - {arg_dict["add_fasta"]}'
                if all([ab_check,id_check]) and not dupl_check:
                    my_database.add_record(acc_nr,sqnc,tax_str,descr)
                    my_database.calculate_content(taxonomy_string=tax_str)
                    my_logger.operation_type = "add_fasta"
                    my_logger._valid_accession = acc_nr
                    my_logger.update_log_dict()
                    my_logger.log_change()
                    my_logger.close_logger()
                else:
                    if (not ab_check) and (not id_check):
                        sys.exit(f'Expected alphabet & Expected id pattern(regex): {my_interface.accession_pattern}|{my_interface.alphabet}')
                    elif not ab_check:
                        sys.exit(f'Expected alphabet: {my_interface.alphabet}')
                    elif not id_check:
                        sys.exit(f'Expected id pattern(regex): {my_interface.accession_pattern}')
                    elif dupl_check:
                        sys.exit(f'Local database already contain a sequence with given accession number.')
                    
        else:
            sys.exit(f'Files with expected type do not exist. Expected:\n {arg_dict["add_fasta"]}.{my_interface.valid_input_type}\n{arg_dict["add_fasta"]}.fasta')

    if arg_dict['add_ncbi']:
        id_check = my_interface.check_format(arg_dict['add_ncbi'])
        dupl_check = my_database.find_id(arg_dict['add_ncbi'])
        if id_check:
            if not dupl_check:
                request_check = my_qry.check_output(arg_dict['add_ncbi'])
                if request_check:
                    output_seq = my_qry.get_fasta(arg_dict['add_ncbi'])
                    sqnc = output_seq[0]
                    acc_nr = arg_dict['add_ncbi']
                    tax_str = "|".join(my_qry.get_taxonomy(arg_dict['add_ncbi']))
                    descr = 'added from ncbi'
                    my_database.add_record(acc_nr,sqnc,tax_str,descr)
                    my_database.calculate_content(taxonomy_string=tax_str)
                    my_logger.operation_type = "add_ncbi"
                    my_logger._valid_accession = acc_nr
                    my_logger.update_log_dict()
                    my_logger.log_change()
                    my_logger.close_logger()
            else:
                sys.exit(f'Local database already contain a sequence with given accession number.')
        else:
            sys.exit(f'Accession id format is not valid. Expected(regex): {my_interface.accession_pattern}')
        
    if arg_dict['add_ncbi_list']:
        file_check = os.path.exists(f'{arg_dict["add_ncbi_list"]}')
        if file_check:
            record_list = list(pd.read_csv(arg_dict['add_ncbi_list'],header=[0]).iloc[:, 0])
            for acc_nr in record_list:
                id_check = my_interface.check_format(acc_nr)
                dupl_check = my_database.find_id(acc_nr)
                if id_check:
                    if not dupl_check:
                        request_check = my_qry.check_output(acc_nr)
                        if request_check:
                            try:
                                output_seq = my_qry.get_fasta(acc_nr)
                                sqnc = output_seq[0]
                                tax_str = "|".join(my_qry.get_taxonomy(acc_nr))
                                descr = 'added from ncbi'
                                my_database.add_record(acc_nr,sqnc,tax_str,descr)
                                my_database.calculate_content(taxonomy_string=tax_str)
                                # my_logger = Logger("add_ncbi",valid_accession=acc_nr)
                                my_logger.operation_type = "add_ncbi"
                                my_logger._valid_accession = acc_nr
                                my_logger.update_log_dict()
                                my_logger.log_change()
                            except ValueError:
                                print(f"{acc_nr} corresponds to sequencing project - no reference sequence attached.")
                                continue
                    else:
                        print(f'Local database already contain a sequence with given accession number: {acc_nr}')
                else:
                    print(f'Accession id format is not valid: {acc_nr} Expected(regex): {my_interface.accession_pattern}')
            my_logger.close_logger()
        else:
            sys.exit(f'File with expected type do not exist. Expected:\n {arg_dict["add_ncbi_list"]}.{my_interface.valid_input_type}')

    if arg_dict['exp_fasta']:
        file_check = os.path.exists(f'{arg_dict["exp_fasta"]}')
        if file_check:
            record_list = list(pd.read_csv(arg_dict['exp_fasta'],header=[0]).iloc[:, 0])
            out_dict = {}
            for acc_nr in record_list:
                exists = my_database.find_id(acc_nr)
                if exists:
                    out_dict[acc_nr] = my_database.export_fasta(acc_nr)
                else:
                    print(f'Sequence with given accession number is not stored in local database: {acc_nr}')
            with open(f'exported_sequences.fasta', 'w') as handle:
                SeqIO.write(out_dict.values(), handle, 'fasta')
            sys.exit("Exported available sequences to exported_sequences.fasta file")
        else:
            sys.exit(f'File with expected type do not exist. Expected:\n {arg_dict["exp_fasta"]}.{my_interface.valid_input_type}')

    if arg_dict['exp_meta']:
        file_check = os.path.exists(f'{arg_dict["exp_meta"]}')
        if file_check:
            record_list = list(pd.read_csv(arg_dict['exp_meta'],header=[0]).iloc[:, 0])
            out_df = pd.DataFrame(columns=["accession_number", "taxonomy_string"])
            for acc_nr in record_list:
                exists = my_database.find_id(acc_nr)
                if exists:
                    record = my_database.export_tax(acc_nr)
                    out_df = out_df.append(record)
                else:
                    print(f'Sequence with given accession number is not stored in local database: {acc_nr}')
            out_df.to_csv(f'exported_taxonomy.csv',header=True, index=False)
            sys.exit("Exported available taxonomies to exported_taxonomy.csv file")
        else:
            sys.exit(f'File with expected type do not exist. Expected:\n {arg_dict["exp_meta"]}.{my_interface.valid_input_type}')

    if arg_dict['exp_records']:
        file_check = os.path.exists(f'{arg_dict["exp_records"]}')
        if file_check:
            record_list = list(pd.read_csv(arg_dict['exp_records'],header=[0]).iloc[:, 0])
            out_df = pd.DataFrame(columns=["accession_number", "taxonomy_string"])
            out_dict = {}
            for acc_nr in record_list:
                exists = my_database.find_id(acc_nr)
                if exists:
                    record = my_database.export_tax(acc_nr)
                    out_df = out_df.append(record)
                    out_dict[acc_nr] = my_database.export_fasta(acc_nr)
                else:
                    print(f'Sequence with given accession number is not stored in local database: {acc_nr}')
            out_df.to_csv(f'exported_taxonomy.csv',header=True, index=False)
            with open(f'exported_sequences.fasta', 'w') as handle:
                SeqIO.write(out_dict.values(), handle, 'fasta')
            sys.exit("Exported available taxonomies to exported_taxonomy.csv file\nExported available sequences to exported_sequences.fasta file")
        else:
            sys.exit(f'File with expected type do not exist. Expected:\n {arg_dict["exp_records"]}.{my_interface.valid_input_type}')

    if arg_dict['rm_record']:
        exists = my_database.find_id(arg_dict['rm_record'])
        if exists:
            tax = my_database.find_tax(arg_dict['rm_record'])
            my_database.rm_record(arg_dict['rm_record'])
            my_database.calculate_content(taxonomy_string=tax, reduce=True, accession=arg_dict['rm_record'])
            my_logger.operation_type = "rm_record"
            my_logger._valid_accession = arg_dict['rm_record']
            my_logger.update_log_dict()
            my_logger.log_change()
            my_logger.close_logger()
        else:
            sys.exit('Sequence with given accession number is not stored in local database.')

    if arg_dict['ch_header']:
        old_id = arg_dict['ch_header'].split(",")[0]
        new_id = arg_dict['ch_header'].split(",")[1]
        id_check_old = my_interface.check_format(old_id)
        id_check_new = my_interface.check_format(new_id)
        id_exists_new = my_database.find_id(new_id)
        exists = my_database.find_id(old_id)
        if id_check_old:
            if exists:
                if id_check_new:
                    if not id_exists_new:
                        tax = my_database.find_tax(old_id)
                        my_database.write_id(old_id, new_id)
                        my_logger.operation_type = "ch_header"
                        my_logger._valid_accession = old_id
                        my_logger._old_accession = old_id
                        my_logger._new_accession = new_id
                        my_database.calculate_content(recalc_accessions=True, accession=old_id, new_accession=new_id, taxonomy_string=tax)
                        my_logger.update_log_dict()
                        my_logger.log_change()
                        my_logger.close_logger()
                    else:
                        sys.exit(f'New id already exists in the database: {new_id}')
                else:
                    sys.exit(f'New id format is not valid: {new_id} Expected(regex): {my_interface.accession_pattern}')
            else:
                sys.exit(f'Sequence with given accession number is not stored in local database: {old_id}')
        else:
            sys.exit(f'Old id format is not valid: {old_id} Expected(regex): {my_interface.accession_pattern}')

    if arg_dict['ch_tax']:
        id = arg_dict['ch_tax'].split(",")[0]
        new_tax = arg_dict['ch_tax'].split(",")[1]
        id_check = my_interface.check_format(id)
        exists = my_database.find_id(id)
        if id_check:
            if exists:
                old_tax = my_database.find_tax(id)
                my_database.write_tax(id,new_tax)
                my_database.calculate_content(old_tax, accession=id, reduce=True)[2]
                my_database.calculate_content(taxonomy_string=new_tax)[2]
                my_logger.operation_type = "ch_tax"
                my_logger._valid_accession = id
                my_logger._old_taxonomy = old_tax
                my_logger._new_taxonomy = new_tax
                my_logger.update_log_dict()
                my_logger.log_change()
                my_logger.close_logger()
            else:
                sys.exit(f'Sequence with given accession number is not stored in local database: {id}')
        else:
            sys.exit(f'Accession format is not valid: {id} Expected(regex): {my_interface.accession_pattern}')

    if arg_dict['view_data']:
        if len(pd.read_csv(my_database.taxonomy_file,header=[0])['accession_number']) > 0:
            count_data = my_database.calculate_content()
            my_plotter = Plotter(count_data[0],count_data[1],count_data[2],count_data[3])
            my_plotter.display()
        else:
            sys.exit(f'Local database does not exist yet or is empty.')
