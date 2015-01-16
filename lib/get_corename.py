# get_corename

import re

def get_corename(protein_id,is_no_corename):
    if is_no_corename:
        return (protein_id)

    m=re.search(r'(.+)(\.[^.]+$)', protein_id)
    if m.group():
        (corename, no_of_isoform) = (m.group(1), m.group(2))  #e.g. AT1G123450.1 => ["AT1G123450", "1"]
        return (corename)
        #(corename, no_of_isoform)=protein_id.split('.')  #e.g. AT1G123450.1 => ["AT1G123450", "1"]

    m=re.search(r'(.+)( [^ ]+$)', protein_id)
    if m.group():
        (corename, others) = (m.group(1), m.group(2))  #e.g. "AT1G123450 123" => ["AT1G123450", "123"]
        return (corename)

    corename = protein_id
    return (corename)

