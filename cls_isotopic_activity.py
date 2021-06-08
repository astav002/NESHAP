import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import json

class cls_isotopic_activity:
    

    def __init__(self) -> None:
        self.define_constants()

    def define_constants(self):

        print()
        print("Initialize isotopic activity parameters")
        print()

        iso_file = open(os.path.join(os.getcwd(), "isotopes.json"), "r")
        self.iso_const = json.load(iso_file)

        # self.iso_const = {
        #     "N-13":{
        #         "t_half_s":600, 
        #         "production":1.1e8 + 4.9e6,
        #         "prod_type": "g"
        #     },
        #     "Ar-41":{
        #         "t_half_s":6.59e3,
        #         "production":0,
        #         "prod_type": "n_thermal"
        #     },
        #     "C-11":{
        #         "t_half_s":1.23e3,
        #         "production":1.1e7+4.5e6,
        #         "prod_type": "g"

        #     },
        #     "H-3":{
        #         "t_half_s":3.88e8,
        #         "production":7.1e6+1.5e7,
        #         "prod_type": "g"
        #     },
        #     "Be-7":{
        #         "t_half_s":4.67e6,
        #         "production":1.1e6+4.5e6,
        #         "prod_type": "g"
        #     },        
        #     "O-15":{
        #         "t_half_s":1.26e2,
        #         "production":5.6e7+4.2e6,
        #         "prod_type": "g"
        #     },           
        #     "Cl-38":{
        #         "t_half_s":2.22e3,
        #         "production":6.8e5,
        #         "prod_type": "g"
        #     },
        #     "Cl-39":{
        #         "t_half_s":3.30e3,
        #         "production":8.5e6,
        #         "prod_type": "g"
        #     }                                  
        # }

        # ventillation rate is also used in cls_air_dat but redefined again here for convenience
        # ventilation is converted to L/s to calculated effected lambda = lambda_decay + lambda_V
        # lambda_V = F / V; we also need surface area and path length for thermal neutron production

        hall_file = open(os.path.join(os.getcwd(), "hall_constants.json"), 'r')
        self.hall_const = json.load(hall_file)

        # self.hall_const = {
        #     "hall_a":{
        #         "volume_L":4.20e4,
        #         "ventilation_cfm":1000,
        #         "ventilation_L_s": 28316.847 / (60 * 1000),
        #         "path_length": 26,
        #         "surface_area": 7576
        #     },
        #     "hall_b":{
        #         "volume_L":1,
        #         "ventilation_cfm":1000,
        #         "ventilation_L_s": 28316.847 / (60 * 1000)
        #     },
        #     "hall_c":{
        #         "volume_L":3.00e4,
        #         "ventilation_cfm":1000,
        #         "ventilation_L_s": 28316.847 / (60 * 1000),
        #         "path_length": 27,
        #         "surface_area": 5925                
        #     },
        #     "hall_d":{
        #         "volume_L":4.20e4,
        #         "ventilation_cfm":1000,
        #         "ventilation_L_s": 28316.847 / (60 * 1000),
        #         "path_length": 26,
        #         "surface_area": 7576                
        #     },
        #     "BSY":{
        #         "volume_L":1.82e4,
        #         "ventilation_cfm":1000,
        #         "ventilation_L_s": 28316.847 / (60 * 1000),
        #         "path_length": 3,
        #         "surface_area": 20000                
        #     }                                                

        # }

        # insert hall decay constant for F/V using the flow rate and location volume
        self.update_w_decay_consts()


    def update_w_decay_consts(self):
        """
        update_w_decay_consts

        Function to add lambda_decay and lambda_volume to the dicts with isotope and hall constants
        there are used in the calculated of lambda_effective for activity calculations
        """
        

        # add the lamba_flow
        for hall in self.hall_const:
            
            self.hall_const[hall]["lambda_v_s"] = (self.hall_const[hall]["ventilation_L_s"] / 
                                                    self.hall_const[hall]["volume_L"])
          

        # now add the lamba_decay
        for isotope in self.iso_const:
            self.iso_const[isotope]["lambda_d_s"] = np.log(2) / self.iso_const[isotope]["t_half_s"]                                                   


    def calculate_activity(self, n13_activity, hall, isotope, t_run_s):
        """
        calculate_activity

        Use the production ratios and isotope constants for a given hall to determine activity
        for a single isotope activity
        """

        # depending on production mechanism, the production value changes in calculation
        # but the overal activity calculation is the same
        if (self.iso_const[isotope]["prod_type"] != "n_thermal"):

            prod_iso = self.iso_const[isotope]["production"]            

        else:

            mult_const = 1.79e5
            vol = self.hall_const[hall]["volume_L"]
            surf = self.hall_const[hall]["surface_area"]
            path = self.hall_const[hall]["path_length"]

            prod_iso = mult_const * vol / (surf * path)


        # calculate the effective decay constant lambda_decay + lambda_exhaust
        lam_iso_eff = self.iso_const[isotope]["lambda_d_s"] + self.hall_const[hall]["lambda_v_s"] 
        lam_n13_eff = self.iso_const["N-13"]["lambda_d_s"] + self.hall_const[hall]["lambda_v_s"]        

        prod_n13 = self.iso_const["N-13"]["production"] 


        lam_n13_d = self.iso_const["N-13"]["lambda_d_s"]


        lam_iso_d = self.iso_const[isotope]["lambda_d_s"]
            

        const_term = prod_iso / prod_n13 * (lam_iso_d / lam_n13_d) * np.power(lam_n13_eff / lam_iso_eff,2 )

        act = (const_term * n13_activity * 
            (t_run_s * lam_iso_eff - (1 - np.exp(-lam_iso_eff * t_run_s))) / 
            (t_run_s * lam_n13_eff - (1 - np.exp(-lam_n13_eff * t_run_s)))
            )                

        return act




if ( __name__ == "__main__"):
    iso = cls_isotopic_activity()


    hall = "hall_a"
    n13_hall_activity = 4.08e5
    t_run_s = 6886800




    for isotope in iso.iso_const:
        if (isotope in iso.iso_const):
            t_half = iso.iso_const[isotope]["t_half_s"]
            prod = iso.iso_const[isotope]["production"]
            iso.iso_const[isotope]["lambda_d_s"]
 

            print("Isotope {}, t_half (seconds): {:2.2e} (lambda_d {:2.2e}), lambda_v {:2.2e},  P: {:2.2e}".format(
                    isotope, 
                    t_half, 
                    iso.iso_const[isotope]["lambda_d_s"],
                    iso.hall_const[hall]["lambda_v_s"],
                    prod))

            print("{}, {} Actvity: {:2.2e} uCi".format(hall, isotope, iso.calculate_activity(n13_hall_activity,hall, isotope, t_run_s)))


