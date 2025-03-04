import numpy as np
import re
import os

class XYParser:
    def __init__(self, file):
        self.pattern = re.compile(r'(\w+): "([^"]+)"')
        self.data_raw_t, self.meta = self.parse(file)

    def get_data(self):
        return self.data_raw_t, self.meta
    
    def parse(self, file):
        meta = {}
        data_raw_t = {}
        
        if isinstance(file, str):
            open_file = open(file, "r")
            ds_name = os.path.basename(file).split(".")[0]
        else:
            open_file = file
            ds_name = "byteio_data"
        
        with open_file as f:
            first_line = f.readline().strip()
            matches = self.pattern.findall(str(first_line))
            for key, value in matches:
                try:
                    value = float(value)
                except ValueError:
                    pass
                meta[key] = value
            
            if isinstance(file, str):
                data_raw_t[ds_name] = np.loadtxt(file, skiprows=1)
            else:
                data_raw_t[ds_name] = np.loadtxt(f, skiprows=1)
        
        return data_raw_t, meta
