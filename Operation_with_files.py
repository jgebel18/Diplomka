import numpy as np
import matplotlib.pyplot as plt
import h5py

import os
#I added this class for shorter time interval to vizualize our results
#from my work i found package h5py and i am using it for saving data accros the iteration process

class Operation_with_files:
    #Like in other classes here
    # i am adding the constructor with some global parameters
    def __init__(self, filename, ):
        self.Files_Path='Files'
        self.filename= os.path.join(self.Files_Path,filename)

#This method is for reading data from
    # datasets with defined name
    def Read_file(self ,nazev):
        with h5py.File(self.filename,'r')as f:
            data= f[(nazev)][:]
        return data

    #This function is for write data into file
    # by dataset with defined name or rewrite this
    # dataset to other data
    def Write_to_file(self, nazev, Data):
        with h5py.File(self.filename, 'a') as f:
            if nazev in f:
                f[nazev][:] = Data  # přepis obsahu
            else:
                f.create_dataset(nazev, data=Data)






