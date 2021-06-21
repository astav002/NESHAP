
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import json

import cls_isotopic_activity as cia

#1 cmf = 28316.847 ml/min
# multiply by the following to convert to ml/min
cfm_to_mlmin = 28316.847

class air_dat():
    
    def __init__(self, **kwargs):
        
        # pull the column mappings from the column maps json file
        col_maps = open(os.path.join(os.getcwd(), "column_maps.json"))
        maps = json.load(col_maps)
        self.monitor = maps["monitor"]
        self.current = maps["current"]
        self.energy = maps["energy"]
        self.current_limit = maps["current_limit"]
        self.air_threshold = maps["air_threshold"]
        self.joins = maps["join"]
        self.current_units = maps["current_unit"]
        
        # self.monitor = {"hall_a":"AirMon_A", 
        #                      "hall_c":"AirMon_C", 
        #                      "hall_b":"AirMon_B", 
        #                      "hall_d":"AirMon_D", 
        #                      "L":"AirMon_L", 
        #                      "SL":"AirMon_SL", 
        #                      "BSY":"AirMon_BSY"}
        
        # self.current = {"hall_a":"C_A(uA)", 
        #                      "hall_c":"C_C(uA)", 
        #                      "hall_b":"C_B(uA)", 
        #                      "hall_d":"C_D(nA)", 
        #                      "BSY":"C_BSY(uA)"}
        
        # self.energy = {"hall_a":"E_A(MeV)", 
        #                "hall_c":"E_C(MeV)", 
        #                "hall_b":"E_B(MeV)", 
        #                "hall_d":"E_C(MeV)", 
        #                "BSY":"E_A(MeV)"}
        
        if ("threshold" in kwargs):
            self.threshold = kwargs["threshold"]
        else:
            self.threshold = 1e-6
            
            
        if ("ventilation" in kwargs):
            self.ventilation = kwargs["ventilation"]
        else:
            self.ventilation = 1000    
            
        print("Ventilation Rate: {} cfm".format(self.ventilation))
            


    def get_data(self, fle_name="./ans2020_Jul29_Sep22_2.airdat.txt", dt_start="07/29/2021 00:00:00"):
        self.fle_name = fle_name
        
        
        air_df = pd.read_csv(fle_name, delimiter="\s+")

        dte_start = datetime.datetime.strptime(dt_start, "%m/%d/%Y %H:%M:%S")
        self.dte_start = dte_start

        
        air_df["DATE_TIME"] = air_df["Time(min)"].apply(lambda x: datetime.timedelta(minutes=x) + dte_start)
        air_df    
        return air_df
    
    def get_date_df(self, air_df):
        """
        set_date_df(self, air_df)
        
        return a dataframe with just the DATE_TIME field for use in merging
        """
        dt_df = pd.DataFrame(air_df["DATE_TIME"])
        return dt_df

    def join_cols(self, air_df):
        # we can join columns into a single extra column for later use
        # example join 'IPM1C01', 'IPM1C02', 'IPM2C01', 'IPM2C02', ... into BSY_CUR_2(uA)

        for col in self.joins:
            print("Populate {} with max of {}".format(col, self.joins[col]))
            air_df[col] = air_df[self.joins[col]].apply(lambda x: max(x), axis=1)

        return air_df





    def apply_threshold(self, air_df, key="hall_c", sigma_method=True, repl_mean=True, replace_zero=True):
        print()
        air_df_sub = air_df[self.monitor[key]]
        #air_df = air_df[air_df_sub < self.threshold]
        threshold = self.air_threshold[key]
        print("Applying threshold to data: {:2.2e} uCi/ml".format(threshold))
        
        print("length before thresholding: {}".format(len(air_df)))
        test_res = air_df_sub > threshold
        
        print("Reject {} using threshold".format(len(test_res[test_res == True])))
        
        
        air_df.loc[test_res, [self.monitor[key]]] = 0
        
        # if we choose we can reject the mean + 3 sigma values once the 
        if (sigma_method):
            print("Replace greater than 3 sigma with mean value")
            threshold = (air_df_sub[air_df_sub <= threshold].std() * 3 + 
                         air_df_sub[air_df_sub <= threshold].mean())
            
            sigma_test_res = air_df_sub > threshold
            print("Found {} values using > 3*sigma".format(len(test_res[sigma_test_res == True])))
            
            air_df.loc[sigma_test_res, [self.monitor[key]]] = 0
            
            test_res = test_res | sigma_test_res
            
        
        
        # Replace our thresholded rejection with the mean of the remainder of the dataset
        if (repl_mean):
            print("Replacing thresholded rejection with mean value")
            mn = air_df[self.monitor[key]].mean()

            print("Mean value replacement: {:2.2e} uCi/ml".format(mn))

            air_df.loc[test_res, [self.monitor[key]]] = mn

            
            
        if (replace_zero):
            
            threshold = ( air_df_sub[air_df_sub <= self.threshold].mean() - 
                         air_df_sub[air_df_sub <= self.threshold].std() * 3)            
            
            print("Replace less than 3 sigma with mean value")
            mn = air_df[self.monitor[key]].mean()
            test_res = air_df[self.monitor[key]] < threshold
            print("Found {} zero values using < 3*sigma".format(len(test_res[test_res == True])))
            air_df.loc[test_res, [self.monitor[key]]] = mn         
        

        print("length after thresholding: {}".format(len(air_df)))

        return air_df
    
    
    def get_loc_bkg(self, air_df, ttle='hist', key="hall_c", pad=0, show_hist=True, current_limit=1, clean_bkg=True):
        # get the background data for the location using a padding following current-on to allow for decay
        print()
        print("Getting background data based on current limit {:2.2e}".format(current_limit))
        
        full_len = len(air_df[self.current[key]])
        
        print("Full Data length: {}".format(full_len))        
        active_cur = air_df[air_df[self.current[key]] > current_limit][self.current[key]]
        
        idx = active_cur.index.values

        print("Original Current-On Data length: {}".format(len(air_df.loc[idx])))   
        
        cur_index = air_df.index.isin(idx)
        
        for i in range(1,pad):
            cur_index = cur_index | air_df.index.isin(idx+i)




        bkg_index = ~cur_index   
                
        idx = np.sort(idx)
        print("Current-on Data length after pad={}: {}".format(pad, len(air_df[cur_index])))

        # current steps are taken from when current is actually on; current-on data set is active_current+pading set
        current_steps = len(active_cur)
        

        
        bkg = air_df.loc[bkg_index, ["DATE_TIME", 
                                   self.current[key], 
                                   self.monitor[key], 
                                   self.energy[key]]]
        
        # set any zeros to average background, helps with tail points that seem to alway be 0"s
        if (clean_bkg):
            bkg.loc[bkg[self.monitor[key]] <=1e-8, [self.monitor[key]]] = bkg[self.monitor[key]].mean()
                
        
        cur = air_df.loc[cur_index, ["DATE_TIME", 
                                   self.current[key], 
                                   self.monitor[key], 
                                   self.energy[key]]]        
        
        
        
        full_data = air_df[["DATE_TIME", 
                                   self.current[key], 
                                   self.monitor[key], 
                                   self.energy[key]]]   
        
        if (show_hist):
            plt.figure()
            fig, ax = plt.subplots()
            bkg[self.monitor[key]].hist(alpha=0.5, label="bkg_"+key, ax=ax, bins=50)
            cur[self.monitor[key]].hist(alpha=0.5, label="monitor_"+key, ax=ax, bins=50)
            ax.set_title("Air Monitor Distribution: {}".format(key))
            ax.set_xlabel("Concentration ($\mu/ml$)")
            ax.ticklabel_format(axis='x', style='scientific', scilimits=(0,0))
            ax.set_ylabel
            ax.legend()  
            plt.savefig(os.path.join(os.getcwd(), ttle + "_" + key+"_bkg_cur_histograms.jpg"))
            plt.close()
        
        
        return full_data, cur, bkg, current_steps

    
    def normalize_get_net(self, dt_df, cur, bkg, key, no_neg=True, plot=True, ttle="plt", 
        use_mn_bkg=False, use_val_bkg=False, usr_mn=0., interpolate=False, split_mean=[]):


        print()
        bkg = bkg.rename(columns={self.monitor[key]:self.monitor[key]+"_Bkg"} )
        
        if (use_mn_bkg):

            bkg = pd.DataFrame(bkg).set_index("DATE_TIME").resample("30Min").asfreq()
            bkg_mn = bkg[self.monitor[key]+"_Bkg"].mean()
            print("Current-on background set to mean {:2.2e} uCi/ml".format(bkg_mn))
            bkg[self.monitor[key]+"_Bkg"] = bkg[self.monitor[key]+"_Bkg"].replace(np.nan,bkg_mn)
            
        elif (use_val_bkg):

            bkg = pd.DataFrame(bkg).set_index("DATE_TIME").resample("30Min").asfreq()
            print("Current-on background set to user value {:2.2e} uCi/ml".format(usr_mn))
            bkg[self.monitor[key]+"_Bkg"] = bkg[self.monitor[key]+"_Bkg"].replace(np.nan,usr_mn)  

        elif ((len(split_mean) > 0 ) and not (use_mn_bkg or use_val_bkg)):
            # split the data set according to the dates within the analysis_configuration.json file
            # then we'll set the background values to the mean within these ranges

            bkg = pd.DataFrame(bkg).set_index("DATE_TIME").resample("30Min").asfreq()
            bkg = bkg.reset_index()
            dt_mn = []

            for dt in split_mean:
                dt_mn.append(datetime.datetime.strptime(dt, "%m/%d/%Y %H:%M:%S"))

            for i in range(len(split_mean)):
                bkg['DATE_TIME'] = pd.to_datetime(bkg['DATE_TIME'])
                if (i==0):
                    print("firts")
                    temp_mn = bkg.loc[bkg['DATE_TIME'] < dt_mn[i], self.monitor[key]+"_Bkg"].mean()
                    bkg.loc[bkg['DATE_TIME'] < dt_mn[i], self.monitor[key]+'_Bkg'] = temp_mn
                    bkg.loc[bkg['DATE_TIME'] < dt_mn[i], self.monitor[key]+'_Bkg'] = (
                        bkg.loc[bkg['DATE_TIME'] < dt_mn[i], self.monitor[key]+'_Bkg'].replace(np.nan,temp_mn))

                elif (i==len(split_mean)):
                    temp_mn = bkg.loc[bkg['DATE_TIME'] >= dt_mn[i], self.monitor[key]+"_Bkg"].mean()

                    bkg.loc[bkg['DATE_TIME'] >= dt_mn[i], self.monitor[key]+"_Bkg"] = temp_mn

                    bkg.loc[bkg['DATE_TIME'] >= dt_mn[i], self.monitor[key]+"_Bkg"] = (
                        bkg.loc[bkg['DATE_TIME'] >= dt_mn[i], self.monitor[key]+"_Bkg"].replace(np.nan,temp_mn))                    
                else:
                    temp_mn = bkg.loc[(bkg['DATE_TIME'] >= dt_mn[i-1]) & (bkg['DATE_TIME'] < dt_mn[i]),
                         self.monitor[key]+"_Bkg"].mean()

                    bkg.loc[(bkg['DATE_TIME'] >= dt_mn[i-1]) & (bkg['DATE_TIME'] < dt_mn[i]),
                         self.monitor[key]+"_Bkg"] = temp_mn       

                    bkg.loc[(bkg['DATE_TIME'] >= dt_mn[i-1]) & (bkg['DATE_TIME'] < dt_mn[i]),
                         self.monitor[key]+"_Bkg"] = (bkg.loc[(bkg['DATE_TIME'] >= dt_mn[i-1]) & (bkg['DATE_TIME'] < dt_mn[i]),
                         self.monitor[key]+"_Bkg"].replace(np.nan,temp_mn))       

            #bkg = pd.DataFrame(bkg).set_index("DATE_TIME").resample("30Min").interpolate(method='linear')
                
            
            
            #print("Current-on background set with split: {:2.2e}, {:2.2e} uCi/ml on date {}".format())

        else:

            bkg = pd.DataFrame(bkg).set_index("DATE_TIME").resample("30Min").interpolate(method="linear")
        
        bkg_df = dt_df.merge(bkg, left_on="DATE_TIME", right_on="DATE_TIME", how="left")
        
        cur = cur.rename(columns={self.monitor[key]:self.monitor[key]+"_Cur"} )
        #cur = pd.DataFrame(cur).set_index("DATE_TIME").resample("30Min").interpolate(method="linear")
        cur_df = dt_df.merge(cur, left_on="DATE_TIME", right_on="DATE_TIME", how="left") 
        
        #dt_df[self.monitor[key]+"_Cur"] = cur_df[self.monitor[key]+"_Cur"]
        
        
        dt_df[self.monitor[key]+"_Bkg"] = bkg_df[self.monitor[key]+"_Bkg"]
        dt_df = dt_df.merge(cur_df, how="left", left_on="DATE_TIME", right_on="DATE_TIME")
        
        
        
        dt_df["net_"+key] = dt_df[self.monitor[key]+"_Cur"] - dt_df[self.monitor[key]+"_Bkg"]
        
        if (no_neg):
            dt_df.loc[dt_df["net_"+key] < 0, ["net_"+key]] = 0
        
        if (plot):
            plt.figure()
            fig, ax = plt.subplots()
            dt_df.plot(x="DATE_TIME", y="net_"+key, ax=ax)
            ax.set_title(key + " Net Concentration")
            plt.xticks(rotation=45)            
            plt.savefig(os.path.join(os.getcwd(), ttle + "_" + key+"_net_concentration.jpg"))
            plt.close()
        
        return dt_df
    
    def total_activity(self, full_df, key, steps, ventilation=1000, t_step_min=60):
        mins = t_step_min
        concentration = full_df["net_"+key].sum()
        mn_concentration = full_df["net_"+key].mean()
        
        # ventilation is specific to the hall/location and passed as argument
        total = concentration * mins * cfm_to_mlmin * ventilation
       

        print("Minutes run: {:2.2e}".format(mins*steps))
        print("Fraction of a year: {:2.2f}".format(mins*steps / (365.24*24*60)))
        print("Sum of concentration: {:2.2e} uCi/ml".format(concentration))
        print("Average Concentration {:2.2e} uCi/ml".format(mn_concentration))
        print("Total: {:2.2e} uCi".format(total))
        
        return total

    def power_calculation(self, df, hall, plot=True, ttle='plt'):
        """
        Simple power calculation using energy and current over the data range
        """
        
        cur_unit = self.current_units[hall]
        cur_hall_key = self.current[hall]
        ene_hall_key = self.energy[hall]
  

        if (cur_unit != 'none'):
            print("Calculating power for {}".format(hall))
            power = (df[cur_hall_key] * df[ene_hall_key]).sum()
            

            if (cur_unit == "nA"):
                pow_unit = "mW-h"
            if (cur_unit == "uA"):
                pow_unit = "W-h"

            print("     Power = {:2.2e} ({})".format(power, pow_unit))
        else:
            return

    
    def normalize(self, df, key, net_set=True, ttle='plt', use_net_res=True):
        net_key = "net_" + key
        cur_hall_key = self.current[key]
        ene_hall_key = self.energy[key]
        air_hall_key = self.monitor[key]
        
        if (net_set):
            air_hall_key = air_hall_key + "_Cur"
        
        norm_c = df[cur_hall_key]/ df[cur_hall_key].max()
        norm_p = (df[ene_hall_key]*df[cur_hall_key])/ (df[ene_hall_key]*df[cur_hall_key]).max()
        
        if (use_net_res):
            norm_air = df[net_key]/ df[net_key].max()
        else:
            norm_air = df[air_hall_key] / df[air_hall_key].max()

        df["norm_"+key] = norm_c
        df["norm_air_"+key] = norm_air

        #full.plot(x="DATE_TIME", )
        fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(20, 10))
        ax[0][0].plot(df["DATE_TIME"], norm_c, "o-", label="Normalized Current", alpha=0.5)
        ax[0][0].plot(df["DATE_TIME"], norm_air, "o-", label="Normalized Concentration", alpha=0.5)
        ax[0][0].set_title(key+ " Normalized Current and Net Concentration")
        ax[0][0].legend()

        ax[0][1].plot(df["DATE_TIME"], norm_p, label="Normalized Power")
        ax[0][1].plot(df["DATE_TIME"], norm_air, label="Normalized Concentration")
        ax[0][1].set_title(key + " Normalized Power and Net Concentration")
        ax[0][1].legend()
        
        ax[1][0].plot(norm_p, norm_air, "x", label="Normalized Airborne")
        ax[1][0].set_title(key+ " Airborned Concentration vs. Current")    
        ax[1][0].legend()        


        ax[1][1].plot(norm_p, norm_air, "x", label="Normalized Airborne")
        ax[1][1].set_title("Airborned Concentration vs. Power")
        ax[1][1].legend()
        plt.savefig(os.path.join(os.getcwd(), ttle + "_" + key+"_current_power_plots.jpg"))
        plt.close()
        
        
        
    
    
    def generate_plots(self, df, key,ttle='gen_plots', net_set=True,):
        print("Generating standard plots")
        
        net_key = "net_" + key
        cur_hall_key = self.current[key]
        ene_hall_key = self.energy[key]
        air_hall_key = self.monitor[key]

        
        if (net_set):
            air_cur_key = air_hall_key + "_Cur"
            air_bkg_key = air_hall_key + "_Bkg"
        
            fig, ax = plt.subplots(figsize=(20,10), nrows=2, sharex=True)
            ax[0].plot(df["DATE_TIME"], df[air_cur_key], "o-", alpha=0.5, label="Current On")
            ax[0].plot(df["DATE_TIME"], df[air_bkg_key], "o-", alpha=0.5, label="Background")
            ax[0].ticklabel_format(style="sci", scilimits=(0,0), axis="y", useOffset=False)
            ax[0].set_title(key + " Calculated Airborne Concentrations")
            ax[0].legend()

            ax[1].plot(df["DATE_TIME"], df[net_key], "o-", alpha=0.5, label="Net Signal")
            ax[1].set_title(key + " Net Airborne Concentration $\mu Ci/ml$")
            ax[1].legend()
            ax[1].ticklabel_format(style="sci", axis="y", scilimits=(0,0), useOffset=False)      
            plt.savefig(os.path.join(os.getcwd(), ttle+ "_" + key+"_bkg_cur_overlay.jpg"))
            plt.close()
            
        else:
            
            fig, ax = plt.subplots(figsize=(20,10))
            ax.plot(df["DATE_TIME"], df[air_hall_key], "o-", alpha=0.5, label="Air Monitor")
            ax.legend()
            ax.set_title("Air Monitor Result")
            plt.savefig(os.path.join(os.getcwd(), ttle + "_" + key+"_air_mon_result.jpg"))
            plt.close()
        
        
        
        
