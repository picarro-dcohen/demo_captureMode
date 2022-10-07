from distutils.ccompiler import gen_lib_options
import pandas as pd
import numpy as np
import h5py

class captureData:
    '''
    from getCaptureData import captureDat

    path = path/to/h5
    cap = captureData(path)
    cap.generate_outputs()
    '''
   


    def __init__(self,path,debug=False):
        self.debug = debug

        #path to private log
        self.df = self.path_to_df(path)
        if self.debug: print(f"loaded h5 file: {self.df.shape}")
        self.feature = 'peak_detector_state'
        self.successful_captures = []
        self.labels = self.get_gas_concs_labels()

    def get_gas_concs_labels(self):
        cols = self.df.columns
        output = []
        for c in cols:
            if c.startswith("broadband_gasConcs"):
                output.append(c)
        return output

    def path_to_df(self, path):
        temp =  h5py.File(path,'r')
        tdf = pd.DataFrame(temp['results'][:])
        return tdf

    def get_first_indicies(self, array):
        #assumes array of indicies with same value i.e. df.loc[df['ValveMask']==3].index
        output = []
        try:
            last = array[0]
        except IndexError:
            print(f"No {self.looking_for} found")
            return
        for i in array:
            if last+1 != i:
                output.append(i)
            last = i
        return output

    def get_last_indices(self, array):
        #assumes array of indicies with same value i.e. df.loc[df['ValveMask']==3].index
        output = []
        try:
            last = array[0]
        except IndexError:
            print(f"No {self.looking_for} found")
            return
        for i in array[1:]:
            if last+1 != i:
                output.append(last)
            last = i
        output.append(i)
        return output

    def get_all_indices(self, target):
        output = self.df[self.df[self.feature]==target].index
        return output

    def get_first_last_indices(self, target):
        array = self.get_all_indices(target)
        firsts = self.get_first_indicies(array)
        lasts = self.get_last_indices(array)
        try:
            return list(zip(firsts,lasts))
        except:
            return
        
    def get_holding_indices(self):
        self.looking_for = "holding"
        output = self.get_first_last_indices(10)
        self.holding = output
        return output

    def get_transitioning_indices(self):
        self.looking_for = "transitioning"
        output = self.get_first_last_indices(9)
        self.transitioning = output

    def pair_transition_holding(self):
        T = self.transitioning
        H = self.holding
        if self.debug:
            try: 
                print(f"  Found {len(T)} transitoning states ")
                print(f"  and {len(H)} holding states")
            except:
                return
        for i,h in enumerate(H):
            temp = {}
            target = h[0]
            target -= 1
            for ii, j in enumerate(T):
                if j[1] == target:
                    temp["transition"] = T[ii]
                    temp ["holding"] = H[i]
                    #TODO: add appropriate timestamp temp[time] = 
                    self.successful_captures.append(temp)
    
    def get_conc_means(self, start, end):
        output = self.df[start:end][self.labels].mean()
        return output
    
    def generate_one_output(self, event):
        tStart, tEnd = event['transition']
        hStart, hEnd = event['holding']

        transition = self.get_conc_means(tStart, tEnd)
        holding = self.get_conc_means(hStart,hEnd)

        tdf = pd.DataFrame(index=self.labels)
        tdf.insert(0, "Transition", transition)
        tdf.insert(1, "Holding", holding)
        tdf.insert(2, "Relative", holding-transition)
        return tdf

    def generate_outputs(self):
        output = []

        self.get_holding_indices()
        self.get_transitioning_indices()
        self.pair_transition_holding()
        for e in self.successful_captures:
            output.append(self.generate_one_output(e))
        if self.debug: print(f" Generated results for {len(output)} capture events")
        return output


if __name__== "__main__":
    paths = ["src/data/NUV1007-20220913-174928Z-DataLog_BTEX_broadband_Private.h5",
            "src/data/NUV1007-20220913-183535Z-DataLog_BTEX_broadband_Private.h5",
            "src/data/NUV1007-20220913-193931Z-DataLog_BTEX_broadband_Private.h5",
            "src/data/NUV1007-20220913-200118Z-DataLog_BTEX_broadband_Private.h5",
            "src/data/NUV1007-20220913-224347Z-DataLog_BTEX_broadband_Private.h5"]

    
    for path in paths:     
        cap = captureData(path, True)
        try:
            cap1 = cap.successful_captures[0]
        except:
            pass
        cap.generate_outputs()


    









