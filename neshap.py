import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
from pandas.plotting import register_matplotlib_converters
import json
import sys


register_matplotlib_converters()


# import class for all NESHAP operations
import cls_air_dat
import cls_isotopic_activity


def main(fle_name, dt_start, halls, run_ttle, res_df, config):
    
    # setup processing parameters from configuration file

    #run_config = open(os.path.join(os.getcwd(), "analysis_configuration.json"), 'r')
    #config = json.load(run_config)

    # convert config into local variables

    # there are fancier ways of doing this, but we'll keep it explicit for interperating
    repl_mean = config["repl_mean"]
    sigma_method = config["sigma_method"]
    replace_zero = config["replace_zero"]
    pad = config["pad"]
    no_neg = config["no_neg"]
    use_mn_bkg = config["use_mn_bkg"]
    use_val_bkg = config["use_val_bkg"]
    usr_mn = config["usr_mn"]
    threshold = config["threshold"]
    split_mean = config["split_mean"]
    interpolate=config["interpolate"]
    smooth_bkg=config["smooth_bkg"]
    smooth_window=config["smooth_window"]

    # repl_mean = True # used in apply_threshold to replace zeros with the overall mean value for all data
    # sigma_method = True# used in apply_threshold to replace  > 3 sigma with the overall mean value for all data
    # replace_zero = False# used in apply_threshold to replace  < 3 sigma with the overall mean value for all data

    # pad = 5 # used in get_loc_bkg to determine how long after beam to associate with current on, in steps not time
    # current_limit = .01 # used in get_loc_bkg to determine if current is on or off

    # # note default method for sbkg correction in normalize_get_net is linear interpolation between single nearest data
    # no_neg = True # used in normalize_get_net to replace negative net values with 0 otherwise keep negatives
    # use_mn_bkg=True # used in normalize_get_net to set background values for current-on to the background mean
    # use_val_bkg=False # used in normalize_get_net to set background values for current-on to a user-specified value
    # usr_mn= 1.87e-7 # used in normalize_get_net when use_val_bkg is true, 0. is the default (no background subtraction)


    # instantiate the class and set threshold

    air = cls_air_dat.air_dat(threshold=threshold)
    air_df = air.get_data(fle_name, dt_start=dt_start)
    air.join_cols(air_df)
    #air_df = air.get_data(fle_name="./ans2020_Jul29_Sep22_2.airdat.txt", dt_start='07/01/2021 00:00:00')


    for hall in halls:

        # pull all of the data loaded from the column_maps.json file to use in function calls etc.
        air_hall_key = air.monitor[hall]
        cur_hall_key = air.current[hall]
        ene_hall_key = air.energy[hall]
        cur_lim_hall_key = air.current_limit[hall]

        # if ( hall == 'hall_d'):
        #     # have to adjust because hall d is specified in nA current, rest are in uA
        #     current_limit = current_limit * 1000




        # output some information for the user regarding the location and dates used
        print("\n****************************************")
        print("\nRunning results for {}".format(hall))
        start =datetime.datetime.strftime(air_df['DATE_TIME'].iloc[0],  '%m/%d/%Y %H:%M:%S')
        end = datetime.datetime.strftime(air_df['DATE_TIME'].iloc[-1],  '%m/%d/%Y %H:%M:%S')
        print("Date Rage {} - {}".format(start, end))

        # get the current (val), background (bkg) data into separate dataframes based on padding and current limits for a hall
        # anything above the current_limit assumes beam on target

        full, val, bkg, current_steps = air.get_loc_bkg(air_df, 
                                                        ttle=run_ttle,
                                                        key=hall, 
                                                        pad=pad, 
                                                        show_hist=True, 
                                                        current_limit=cur_lim_hall_key)           


        # apply theshold to the full dataset, threshold here is air_threshold from column_maps.json

        bkg = air.apply_threshold(bkg, key=hall,
                                    repl_mean=repl_mean,
                                    sigma_method=sigma_method,
                                    replace_zero=replace_zero)

        val = air.apply_threshold(val, key=hall,
                                    repl_mean=repl_mean,
                                    sigma_method=sigma_method,
                                    replace_zero=replace_zero)                                    





        # we setup a blank dataframe with the date range we want to populate with bkg/air monitoring after processing
        dt_df = pd.DataFrame(air_df['DATE_TIME'])


        full = air.normalize_get_net(dt_df, val, bkg,
                                    key=hall, 
                                    no_neg=no_neg, 
                                    ttle=run_ttle,
                                    plot=True, 
                                    use_mn_bkg=use_mn_bkg, 
                                    use_val_bkg=use_val_bkg, 
                                    usr_mn=usr_mn,
                                    interpolate=interpolate,
                                    split_mean=split_mean,
                                    smooth_bkg=smooth_bkg,
                                    smooth_window=smooth_window
                                    )

        t_step_min=60
        
        iso = cls_isotopic_activity.cls_isotopic_activity()
        ventilation = iso.hall_const[hall]["ventilation_cfm"]

        total_uCi = air.total_activity(full, hall, current_steps, 
                                        ventilation=ventilation, t_step_min=t_step_min )

        t_run_s = current_steps * t_step_min * 60

        print("************************")
        print()
        print("Calculating isotopic activities from total")
        print("Running time: {:2.2f} s".format(t_run_s))



        for isotope in iso.iso_const:
            if (isotope in iso.iso_const):
                t_half = iso.iso_const[isotope]['t_half_s']
                prod = iso.iso_const[isotope]['production']
                lambda_d_s = iso.iso_const[isotope]['lambda_d_s']
                lambda_v_s = iso.hall_const[hall]['lambda_v_s']


                iso_act = iso.calculate_activity(total_uCi, hall, isotope, t_run_s)
                res_df.loc[len(res_df.index)] = [run_ttle, dt_start, hall, isotope, iso_act, lambda_d_s, lambda_v_s, prod]
    
                
                print("Isotope {}, t_half (seconds): {:2.2e} (lambda_d {:2.2e}), lambda_v {:2.2e},  P: {:2.2e}".format(
                        isotope, 
                        t_half, 
                        lambda_d_s,
                        lambda_v_s,
                        prod))
                
                print("{} Actvitivity: {:2.2e} uCi".format(isotope, iso_act))       

        #res_df.to_csv(os.path.join(os.getcwd(), run_ttle + "_" + hall+ ".csv"))


        # Generate standard plots
        air.generate_plots(full, hall, ttle=run_ttle, net_set=True)
        air.power_calculation(full, hall,plot=True, ttle=run_ttle)
        print()

    return res_df


if (__name__ == "__main__"):


    if (len(sys.argv) > 1):
        print("\n************************************************")
        print("\nReading in user specified configuration file: {}".format(sys.argv[1]))
        
        run_config = open(sys.argv[1])
        config = json.load(run_config)
        fle_name = config["data_files"]
        dt_start = config["data_starts"]
        halls = config["halls"]
        run_ttle = config["run_titles"]
        output_title = config["output_title"]

    else:
        run_config = open(os.path.join(os.getcwd(), "configurations", "analysis_configuration.json"), 'r')
        config = json.load(run_config)
        fle_name = config["data_files"]
        dt_start = config["data_starts"]
        halls = config["halls"]
        run_ttle = config["run_titles"]
        output_title = config["output_title"]


        # fle_name = ["./ans2020_Jan01_Mar31_wBSY.airdat.txt", "./ans2020_Jul29_Sep22_wBSY.airdat.txt"]
        # dt_start = ['01/01/2020 00:00:00', '07/29/2020 00:00:00']


    # this is the holder for the results
    res_df = pd.DataFrame(columns=["Run", "Start Date", "Hall", "Isotope", "Activity(uCi)", "Lambda_D(s)", "Lambda_V(s)", "P"])

    # usr_fle = input("Enter the file name using relative or full path: \n")

    # if (usr_fle == ""):
    #     fle_name = fle_name
    # else:
    #     fle_name = usr_fle



    for i in range(len(fle_name)):
        res_df = main(fle_name[i], dt_start[i], halls, run_ttle[i], res_df,config)
        
    res_df.to_csv(os.path.join(os.getcwd(), output_title + ".csv"))
