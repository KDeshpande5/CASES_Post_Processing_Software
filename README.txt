-------------------------------------------------------------
CASES_Post_Processing_Software version v1 (19 Jan 2014):
-------------------------------------------------------------

This software has been developed by Kshitija Deshpande at Virginia Tech in June 2012 to post process scintillation data from Connected Autonomous Space Environment Sensor (CASES) Global Positioning System (GPS) software-defined receivers by ASTRA (http://cases.astraspace.net). More details on the processing are in a 2012 Radio Science (2012RS005061) journal paper (citation below). Please cite this paper if you are using this CASES Post Processing software.

K. Deshpande, G. S. Bust, C. R. Clauer, H. Kim, J. E. Macon, T. E. Humphreys, J. A. Bhatti, S. B. Musko, G. Crowley, and A. T. Weatherwax (2012), Initial GPS Scintillation results from CASES receiver at South Pole, Antarctica, Radio Sci., 47, RS5009, doi:10.1029/2012RS005061.

Each code in this software is provided as is. The codes are available for distribution. Please send comments, questions and bug reports to Kshitija Deshpande at kshitija@vt.edu.

The post processing involves high rate and low rate data processing.

1. High rate data processing
Run_Read_Plot_HighRateCASESdata_routines.m 
Calls routines: 
Fn_ReadHighRate_CASESdata.m to read CASES high rate data and
Fn_Plot_HighRate_CASESdata.m to plot CASES high rate data in several ways (all scintillating events, S4 and sigma_phi etc.). 

Needs to input the correct path name and folder name for the CASES data (where the daily .log files are). 

See Fn_ReadHighRate_CASESdata.m for details on signal type. Generally, it is 0 or 2 for AAL-PIP CASES. Thus, i=1 or 3. Also, different plot types that can be created with the plotting routines are mentioned in the same code. Select what is necessary, or just plot everything.


2. Low rate data processing
All .m files need to input the correct path name and folder name for the CASES data (where the daily .log files are). 

a. LowRateCASES_S4Sigphi.m -> plots sigma_phi and S4 for all available low rate data
b. LowRateCASES_StackPlotting.m -> creates 2 stacked plots each for S4 and sigma_phi containing all PRNs.
c. LowRateCASES_TEC_plot.m -> creates TEC plots from available low rate data


The rest of the .m files (AJbutter.m, butterworth_discrete.m, err_chk.m, gps2utc.m, slidingS4SigmaPHI.m, utc2leap.m) are functions required by these codes.


This folder also includes 2 folders: a) CASES_DeflateSoftware/ -- to generate executables (binflate and sbcclient) for deflating the CASES GPS receiver data, and b) Shell_scripts/ -- for shell scripts to help deflate the data. Please refer to the individual read me files inside these directories for further details.

