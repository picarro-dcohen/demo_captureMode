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
        self.path = path
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
        try:
            self.pubChemData = pd.read_csv('src/CID_list.csv')
        except:
            pass 
        return output
    
    def make_multiplier_mask(self):
        multiplier = []
        unit = []
        cid_list = [int(cid.split('_')[-1]) for cid in self.labels]
        for cid in cid_list:
            if cid == 280 or cid == 297:
                multiplier.append(1e6)
                unit.append('ppm')
            elif cid == 962:
                multiplier.append(1)
                unit.append('-')
            else:
                multiplier.append(1e9)
                unit.append('ppb')
            
        return multiplier, unit, cid_list


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
            print(f"  No {self.looking_for} found")
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
            #if self.debug: print(f"Did not find any {self.looking_for}")
            return
        
    def get_holding_indices(self):
        '''
        gets the first and last indicies of all holding indicies
        holding is defined where:
            self.feature = peak_detector_state == 10
        '''
        self.looking_for = "holding"
        output = self.get_first_last_indices(10)
        self.holding = output
        return output

    def get_transitioning_indices(self):
        '''
        gets the first and last indicies of all transitioning indicies
        transitioning is defined where:
            self.feature = peak_detector_state == 9
        '''
        self.looking_for = "transitioning"
        output = self.get_first_last_indices(9)
        self.transitioning = output

    def get_triggered_indices(self):
        '''
        gets the first and last indicies of all triggered indicies
        triggered is defined where:
            self.feature = peak_detector_state == 3
        '''
        self.looking_for = "triggered"
        output = self.get_first_last_indices(3)
        self.triggered = output

    def pair_transition_holding(self):
        trig = self.triggered
        T = self.transitioning #list of first/last indicides 
        H = self.holding
        if self.debug:
            try: 
                print(f"  Found {len(trig)} triggered states ")
                print(f"  Found {len(T)} transitoning states ")
                print(f"  Found {len(H)} holding states")
            except:
                print("Found No transition or holding pairs.")
                return
        if H is None or T is None:
            return

        for i,h in enumerate(H):
            temp = {}
            target = h[0]
            target -= 1
            for ii, j in enumerate(T):
                if j[1] == target:
                    temp["triggered"] = trig[i]
                    temp["transition"] = T[ii]
                    temp ["holding"] = H[i]
                    #TODO: add appropriate timestamp temp[time] = 
                    self.successful_captures.append(temp)

    def pair_triggered_transition(self):
        pass


    
    def get_conc_means(self, start, end):
        output = self.df[start:end][self.labels].mean()
        return output
    
    def generate_one_output(self, event):
        trigStart,trigEnd =event['triggered']
        tStart, tEnd = event['transition']
        hStart, hEnd = event['holding']
        triggered = self.get_conc_means(trigStart,tEnd)
        transition = self.get_conc_means(tStart, tEnd)
        holding = self.get_conc_means(hStart,hEnd)
        multiplier, units, cid_list = self.make_multiplier_mask()

        tdf = pd.DataFrame(index=self.labels)
        tdf.insert(0, "Name", list(map(self.parse_name,self.labels)))
        tdf.insert(1, "Triggered", triggered*multiplier)
        tdf.insert(2, "Transition", transition*multiplier)
        tdf.insert(3, "Holding", holding*multiplier)
        relative = (holding-triggered)*multiplier
        tdf.insert(4, "Relative", relative)
        print(np.sum(triggered),triggered[4], triggered[-1]/np.sum(triggered))
        tdf.insert(5, "Units", units)
        tdf = self.water_to_pct(tdf, triggered, transition, holding, relative)

        tdf = tdf.round(decimals=4)
        return tdf

    def water_to_pct(self,tdf, triggered, transition, holding, relative):
        old_row = tdf.loc['broadband_gasConcs_962']
        tdf.at['broadband_gasConcs_962','Triggered'] = old_row['Triggered']/np.sum(triggered)
        tdf.at['broadband_gasConcs_962','Transition'] = old_row['Transition']/np.sum(transition)
        tdf.at['broadband_gasConcs_962','Holding'] = old_row['Holding']/np.sum(holding)
        tdf.at['broadband_gasConcs_962','Relative'] = old_row['Relative']/np.sum(relative)
        tdf.at['broadband_gasConcs_962','Units'] = 'percent' 
        return tdf
        



    def generate_one_meta(self, event):
        indicies = ['triggered','transition','holding']
        columns = ['start_time','end_time', 'duration','start_index','end_index','N']
        tdf = pd.DataFrame(index=indicies, columns=columns)
        tdf.loc['triggered']['start_time'] = self.df.iloc[event['triggered'][0]]['time']
        tdf.loc['triggered']['end_time'] = self.df.iloc[event['triggered'][1]]['time']
        tdf.loc['triggered']['duration'] = tdf.loc['triggered']['end_time'] - tdf.loc['triggered']['start_time']
        tdf.loc['triggered']['start_index'] = event['triggered'][0]
        tdf.loc['triggered']['end_index'] = event['triggered'][1]
        tdf.loc['triggered']['N'] = event['triggered'][1] - event['triggered'][0]

        tdf.loc['transition']['start_time'] = self.df.iloc[event['transition'][0]]['time']
        tdf.loc['transition']['end_time'] = self.df.iloc[event['transition'][1]]['time']
        tdf.loc['transition']['duration'] = tdf.loc['transition']['end_time'] - tdf.loc['transition']['start_time']
        tdf.loc['transition']['start_index'] = event['transition'][0]
        tdf.loc['transition']['end_index'] = event['transition'][1]
        tdf.loc['transition']['N'] = event['transition'][1] - event['transition'][0]
        
        tdf.loc['holding']['start_time'] = self.df.iloc[event['holding'][0]]['time']
        tdf.loc['holding']['end_time'] = self.df.iloc[event['holding'][1]]['time']
        tdf.loc['holding']['duration'] = tdf.loc['holding']['end_time'] - tdf.loc['holding']['start_time']
        tdf.loc['holding']['start_index'] = event['holding'][0]
        tdf.loc['holding']['end_index'] = event['holding'][1]
        tdf.loc['holding']['N'] = event['holding'][1] - event['holding'][0]
        return tdf


    def get_CID_name(self, cid):
        cid = int(cid)
        temp = self.pubChemData.loc[self.pubChemData['cid'] == cid]['cmpdname'].values[0]
        return temp
  

    def parse_name(self, label):
        cid = int(label.split('_')[-1])
        return self.get_CID_name(cid)


    def generate_outputs(self):
        '''
        generates all outputs from a .dat file.
        1. gets all transition indicies start/stop
        2. '' holding ''
        3. pairs each transition to holding -> a pair is a successful capture
        4. generates simple output table and metadata
        '''
        
        meta = []
        output = []
        self.get_triggered_indices()
        self.get_holding_indices()
        self.get_transitioning_indices()
        self.pair_transition_holding()
        for e in self.successful_captures:
            output.append(self.generate_one_output(e))
            meta.append(self.generate_one_meta(e))
        if self.debug: print(f" Generated results for {len(output)} capture events")
        return meta, output

    def save_outputs(self, data):
        out_prefix = './output/'
        out_suffix = self.path.split('/')[-1].split('.')[0]
        for i, e in enumerate(data):
            e.to_csv(out_prefix+out_suffix+"_"+str(i) )


        



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
        meta, data = cap.generate_outputs()
        cap.save_outputs(data)


    









