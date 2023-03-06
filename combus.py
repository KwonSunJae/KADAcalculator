import numpy as np
import csv
import os


def convert_cell_value(val):
    # single value only. no arrays
    if val.replace('.','').isdigit() or val.replace('-','').isdigit():
        return float(val[0])
    else:
        return str(val[0])


class CommunicationBus:
    def __init__(self, template_file):
        path = os.path.abspath(template_file)
        with open(path, mode ='r') as csv_file:
            reader = csv.reader(csv_file, delimiter=',')
            for row in reader:
                self.__dict__[row[0]] = convert_cell_value(row[1])
    
    def __repr__(self):
        out = ''
        for attr, value in self.__dict__.items():
            out += '{:<25} = {:<12}\n'.format(attr,value)
        return out