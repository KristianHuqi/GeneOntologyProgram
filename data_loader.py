
#PARSING THE goa_human.gaf file
import pandas as pd

class GAF_parser:
    
    def __init__(self,filepath):
        self._filepath = filepath
    
    def load_gaf(self):
        print(f'reading the file from: {self._filepath}')
    
        
    
        columns = ['DB','DB object id', 'DB object symbol', 'Relation','GO ID','Reference','Evidence code',
               'With (or) From','Aspect','DB object name','DB object synonym','DB object type',
               'Taxon','Date','Assigned by','Annotation Extension','Gene product form id']
    
    
        
        data_frame = pd.read_csv(self._filepath, sep = '\t', comment = '!', names = columns, dtype = str)
        
    
        #closing step
        print('reading complete')
        return data_frame
'''
#testing the program:
filepath = '/Users/giovannimunari/Desktop/Bioproject/data/goa_human.gaf' #to open the program
data = GAF_parser(filepath) #our dataset obtained as product of the above function
data_set = data.load_gaf()
#print first 5 and last 5 rows as test
try:
    first_5 = data_set.head()
    last_5 = data_set.tail()
    print(first_5)
    print()
    print(last_5)

except FileNotFoundError:
    print('file non trovato, controlla di aver messo il percorso adeguato oppure che i file siano nel formato giusto')
    
print()
print()
print('--------------------------------------------------------------------------------------')
'''

#PARSING THE go-basic.obo file

#here we can't use pandas as it is not a table so we need to write a manual parser
#what we know is that each index is introduced by " [term] " in the file so we use this to build our parser as a dictionary as a list of lists and then convert it into an actual dataframe
  

class OBO_parser:
    
    def __init__(self, filepath):
        self._filepath = filepath

    def parse_obo(self):
        with open(self._filepath, 'r') as f:
            content = f.readlines()

        dictionary_for_dataset = {}
        term_counter = 0 
        current_key = None
        is_inside_term = False

        
        for line in content:
            line = line.strip()
    
            if line == '':
                continue
            
            if line == '[Term]':
                term_counter += 1
                current_key = f'Term {term_counter}'
                dictionary_for_dataset[current_key] = []
                is_inside_term = True
                continue
    
            
            elif line.startswith('[') and line != '[Term]':
                is_inside_term = False #we stop saving
                continue
    
            
            if is_inside_term and current_key is not None:
                parts = line.split(':', 1) 
                
                if len(parts) == 2:
                    parts[0] = parts[0].strip()
                    parts[1] = parts[1].strip()
                    
                    dictionary_for_dataset[current_key].append([parts[0],parts[1]])
    
        #NOW WE CONVERT THE DICTIONARY TO AN ACTUAL PANDAS DATAFRAME

        print('elaborating data and converting in dataframe...')
        terms_list = []
        rels_list = []

        for term, attributes in dictionary_for_dataset.items():
            term_data = {'go_id': None, 'name': None , 'def': None, 'namespace': None, 'exact_synonym': None,'synonyms': [], 'is_obsolete': False, 'replaced_by': None}
    
            
            for key, value in attributes:
                if key == 'id':
                    term_data['go_id'] = value
                elif key == 'name':
                    term_data[key] = value
                elif key == 'namespace':
                    term_data[key] = value
                elif key == 'exact_synonym':
                    term_data[key] = value
                elif key == 'def':
                    if '"' in value:
                        term_data[key] = value.split('"')[1]
                    else:
                        term_data[key] = value
                        
                elif key == 'synonym':
                    if '"' in value:
                        clean = value.split('"')[1]
                    else:
                        clean = value
                    term_data['synonyms'].append(clean)
                    
                elif key == 'replaced_by':
                    term_data[key] = value
                elif key == 'is_obsolete' and value == 'true':
                    term_data[key] = True
    
        
            for key, value in attributes:
                if key == 'is_a':
                    parent_id = value.split('!')[0].strip()
                    rels_list.append({'child_id': term_data['go_id'], 'parent_id': parent_id, 'relation_type': 'is_a'})
    
            terms_list.append(term_data)
        
        terms_df = pd.DataFrame(terms_list)
        rels_df = pd.DataFrame(rels_list)
    
        
        return terms_df, rels_df 