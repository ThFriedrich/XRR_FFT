

from src.brml_parser import BrmlParser
from src.xy_parser import XYParser

def parse_data(file):
    if file.endswith('.brml'):
        data_raw_t, meta = BrmlParser(file).get_data()
    elif file.endswith('.xy'):
         data_raw_t, meta = XYParser(file).get_data()
    else:
        meta = {}
        data_raw_t = {}
    return data_raw_t, meta

def parse_data_bytesIO(file, type):
    if type == 'brml':
        data_raw_t, meta = BrmlParser(file).get_data()
    elif type == 'xy':
        data_raw_t, meta = XYParser(file).get_data()
    else:
        meta = {}
        data_raw_t = {}
    return data_raw_t, meta