%% Run_Read_Plot_HighRateCASESdata_routines.m 
% This code is a part of CASES_Post_Processing_Software v1
% Copyright (C) 2012 - 2014, Kshitija Deshpande, Virginia Tech. 

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% -----------------------------------------------------------------------

% Run_Read_Plot_HighRateCASESdata_routines.m calls routines: 
% 1. Fn_ReadHighRate_CASESdata.m to read CASES high rate data
% and Fn_Plot_HighRate_CASESdata.m to plot CASES high rate data in several
% different ways (all scintillation events, S4 and sigma_phi etc). 
% 2. LowRateCASES_S4Sigphi_elvfiltered.m, 
% LowRateCASES_StackPlotting_elvfiltered.m and LowRateCASES_TEC_plot.m to
% plot low rate S4, sigma_phi for each PRN, as stack plots and low rate 
% TEC plots.

% % Inputs: 
% 1. Needs to change "path" to enter the correct path of the CASES
% folder and folder name where all the CASES log files are stored.
% folder_path is the complete folder path including the folder name to be
% input to the read high rate data MATLAB function
% (Fn_ReadHighRate_CASESdata) and plotting high rate data function
% (Fn_Plot_HighRate_CASESdata). 
% 2. signal type:
% See Fn_ReadHighRate_CASESdata.m for details on signal type.
% Generally, it is 0 or 2 for AAL-PIP CASES. Thus, i=1 or 3. 
% 3. types of plots
% A = all & scintillating events, all segments
% of data -- processed
% B = S4 and Sigma_phi plots
% C = plot common scintillating times for PRNs
% 
% % Outputs: 
% Different types of plots depending on the user preferences.


close all
clear all

sep = filesep;

%User defined path to the log files of the data.
%Change it for your system
path = '/Users/kshitija/work/GPS_data/'; 

%Change folder name according to the date of CASES data, folder with all
%the log files
folder_name = 'Receiver241/2013/337';
folder_path = [path,folder_name,sep];

Receivername = 'Receiver241';

disp(folder_path)

for i = 1:1:1; % for signal type,
    %see Fn_ReadHighRate_CASESdata.m for details on signal type
    signal_type = i-1;
    Fn_ReadHighRate_CASESdata(folder_path,signal_type);
    for ii = 1:1:3; %Type of plots
        %specify types of plots
        % ii = 1 -> A = all & scintillating events, all segments
        % of data -- processed
        % ii = 2 -> B = S4 and Sigma_phi plots
        % ii = 3 -> C = plot common scintillating times for PRNs
        set_plotCT = ii-1;
        switch set_plotCT
            case 0
                set_plot = 'A';
            case 1
                set_plot = 'B';
            case 2
                set_plot = 'C';
            otherwise
                error('Unknown set.')
        end
        Fn_Plot_HighRate_CASESdata(folder_path,signal_type,set_plot)
    end %set_plot
    %Plot low rate data
        LowRateCASES_S4Sigphi_elvfiltered(folder_path,signal_type,...
            Receivername);
        LowRateCASES_StackPlotting_elvfiltered(folder_path,signal_type,...
            Receivername);
        LowRateCASES_TEC_plot(folder_path,signal_type);
end %signal type
