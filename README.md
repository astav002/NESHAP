# NESHAP

## Analysis configuration
    repl_mean = True # used in apply_threshold to replace zeros with the overall mean value for all data
    sigma_method = True# used in apply_threshold to replace  > 3 sigma with the overall mean value for all data
    replace_zero = False# used in apply_threshold to replace  < 3 sigma with the overall mean value for all data

    pad = 5 # used in get_loc_bkg to determine how long after beam to associate with current on, in steps not time
    current_limit = .01 # used in get_loc_bkg to determine if current is on or off

    # note default method for sbkg correction in normalize_get_net is linear interpolation between single nearest data
    no_neg = True # used in normalize_get_net to replace negative net values with 0 otherwise keep negatives
    use_mn_bkg=True # used in normalize_get_net to set background values for current-on to the background mean
    use_val_bkg=False # used in normalize_get_net to set background values for current-on to a user-specified value
    usr_mn= 1.87e-7 # used in normalize_get_net when use_val_bkg is true, 0. is the default (no background subtraction)

## Constants (Halls, Isotopes, Column Mappings)

### Isotopes
Isotope Constants are stored in a json file "isotope_constants.json".  The derivatio of the production values are shown below, but only the result is stored in the json file :

    iso_const = {
            "N-13":{
                "t_half_s":600, 
                "production":1.1e8 + 4.9e6,
                "prod_type": "g"
            },
            "Ar-41":{
                "t_half_s":6.59e3,
                "production":0,
                "prod_type": "n_thermal"
            },
            "C-11":{
                "t_half_s":1.23e3,
                "production":1.1e7+4.5e6,
                "prod_type": "g"
            },
            "H-3":{
                "t_half_s":3.88e8,
                "production":7.1e6+1.5e7,
                "prod_type": "g"
            },
            "Be-7":{
                "t_half_s":4.67e6,
                "production":1.1e6+4.5e6,
                "prod_type": "g"
            },        
            "O-15":{
                "t_half_s":1.26e2,
                "production":5.6e7+4.2e6,
                "prod_type": "g"
            },           
            "Cl-38":{
                "t_half_s":2.22e3,
                "production":6.8e5,
                "prod_type": "g"
            },
            "Cl-39":{
                "t_half_s":3.30e3,
                "production":8.5e6,
                "prod_type": "g"
            }                                  
        }
