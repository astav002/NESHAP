# NESHAP

Python scripts to automate the processing of JLab data for NESHAPS reporting.  Equations are documented in various TBDs and full SQA is documented within JLab's Docushare system with the Radiation Control Department files.  The bulk of the automation is within the background correction process for determining the N-13 activity from the measured airborne radioactivity levels.  The processing parameters are all controled by an analysis_configuration.json file which can be specified at run time as a command line argument when running the code:

```
python neshap.py ./configuration_files/analysis_configuration.json >> ./output_data/result.txt
```

Requires (or higher version):

* Python 3.6.8 
* Numpy 1.16.4
* Pandas 0.25.1
* Matplotlib 3.1.1



    
