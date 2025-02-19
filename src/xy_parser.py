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
        with open(file, "r") as f:
            first_line = f.readline().strip()
            matches = self.pattern.findall(first_line)
            for key, value in matches:
                try:
                    value = float(value)
                except ValueError:
                    pass
                meta[key] = value
        ds_name = os.path.basename(file).split(".")[0]
        data_raw_t[ds_name] = np.loadtxt(file, skiprows=1)
        return data_raw_t, meta
