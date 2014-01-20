%% LowRateCASES_TEC_plot.m 
% This code is a part of CASES_Post_Processing_Software v1
% Copyright (C) 2012 - 2014, Kshitija Deshpande, Virginia Tech. 

% CASES_Post_Processing_Software v1 is free software: you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% CASES_Post_Processing_Software v1 is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% ----------------------------------------------------------------------

% %% Low rate CASES TEC plotting  


function LowRateCASES_TEC_plot(folder_name,signal_type)

tic


%user defined path to the log files of the data.
folder = dir(folder_name);
s={folder.name}; % read all folders and files in the directory

%File separator for the selected operating system
sep = filesep;

%Specify plot type: plot_type = 1 for eps, anything else for png
plot_type = 0;

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)

%%specify signal type
%signal_type = 0; % 0 = L1C, 1 = L2C etc. 
%               0       GPS_L1_CA       
%               1       GPS_L2_CM       // GPS L2 civil M code
%               2       GPS_L2_CL       
%               3       GPS_L2_CLM       
%               4       GPS_L5_I        
%               5       GPS_L5_Q        
%               6       GPS_L5_IQ       // GPS L5 civil I+Q combined tracking
%               7       GPS_L1_CA_ALT1  // GPS L1 C/A for alternative L1C bank
%               8       CDMA_UHF_PILOT  // Cellular CDMA pilot, I+Q signal
%               9       CDMA_UHF_SYNC   // Cellular CDMA pilot+sync, I+Q signal
switch signal_type
   case 0
      signal = 'L1CA';  %// GPS L1 legacy civil code 
   case 1
      signal = 'L2CM';  %// GPS L2 civil L code
   case 2
      signal = 'L2CL';  %// GPS L2 M+L combined tracking
   case 3
      signal = 'L2CLM';  %// GPS L5 civil in-phase
   case 4
      signal = 'L5I';   %// GPS L5 civil quadrature
   case 5
      signal = 'L5Q';
   case 6
      signal = 'L5IQ';
   case 7
      signal = 'L1CA-ALT1';
   case 8
      signal = 'CDMA-UHF-PILOT';
   case 9
      signal = 'CDMA-UHF-SYNC';
   otherwise
      error('Unknown signal.')
end


%Plotting intialization
fontsz = 10;
marksz = 4;
linestyle= '.';
fontsz_lg = 12; %font size of legends

command = strcat('mkdir',{' '},folder_name,'plots_',signal);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_name,'plots_',signal,sep,'LowRate',sep);
system(cell2mat(command));
command = strcat('mkdir',{' '},folder_name,'prn_files_',signal);
system(cell2mat(command));
outdir = strcat(folder_name,'plots_',signal,sep,'LowRate',sep);


PRN = zeros(1);
ORTW = zeros(1); %ORT weeks
ORTS = zeros(1); %ORT seconds
TEC = zeros(1);
DTEC = zeros(1);

%Specify Input file names
ionfilename = strcat(folder_name,'iono.log');

%% Reading scintillating and reference channels from Ion log file
ionfile = load(ionfilename);
%no dependence on signal type
iondata = ionfile;
%Corresponding time stamps
ort_ion =(iondata(:,2)+iondata(:,3));

% %Find and remove repeated entries
%1st derivative to find the time slots (wherever there is a huge
%difference, that's the BUG. There are repeated entries inside scin.log
ort_scd = diff(ort_ion);
scdn_ind=find(ort_scd<-20);
for ii = 1:1:length(scdn_ind),
    indii = find(ort_ion == ort_ion(scdn_ind(ii)));
    ort_ion(indii(1):indii(end)-1)=0; 
    %This will work well only if there are only 2 repetitions max. ort_ion =
    %ort_ion(ort_ion~=0) will remove the repeated values. But better leave
    %the zeros as they are and recreate the original "iondata" array based
    %on those zeros.
end

%Non-zero time stamps means without any repetition.
nzero_ort = find(ort_ion~=0);
iondata = iondata(nzero_ort,:);

for kk = 1:1:32,
    sv = find(iondata(:,6)==kk);
    data = iondata(sv,:);
    if(isempty(data)==0),
        ortw = data(:,1); % ort week

        %check for valid files by looking at week number
        %GPS week number cannot be greater than 2000
        %If not valid data, come out
        ind_valid = find(ortw > 2000);
        if(length(ind_valid) > 1),
            string = strcat(s(ii),' file has invalid entries for GPS time.');
            disp(string)
            break
        end

        orts = data(:,2) + data(:,3); %ort seconds + fractional seconds

        tec = data(:,4); % TEC
        dtec = data(:,5); % dTEC
        
        prn = data(:,6); %prn

        ORTW = [ORTW, ortw'];
        ORTS = [ORTS, orts'];

        TEC = [TEC, tec'];
        DTEC = [DTEC, dtec'];
        
        PRN = [PRN, prn'];

    end

end %for loop for prns

%Get rid of the leading zero
ORTW = ORTW(2:end);
ORTS = ORTS(2:end);
TEC = TEC(2:end);
DTEC = DTEC(2:end);
PRN = PRN(2:end);

DATA = [ORTW; ORTS; TEC; DTEC; PRN]';
% % struct('Data',DATA);
outfilename = strcat(folder_name,'prn_files_',signal,sep,'Daily_iondata.mat');
save(outfilename,'DATA');


ORTS_BEG = ORTS(1);
ORTS_END = ORTS(end);
L = length(DATA);
ORTS_min = min(ORTS);

% %gps to utc conversion
% [UTC_time, leap_sec, day_of_year]=gps2utc([ORTW(1) ORTS_min], 0);
% %UTC_time = year month day hour minutes seconds
% UTC_init = UTC_time(1,:);
% YEAR = UTC_init(1);
% MONTH = UTC_init(2);
% DAY = UTC_init(3);
% HOUR = UTC_init(4);

% %to completely randomize the colors of stack plot
% RandStream.setDefaultStream ...
%      (RandStream('mt19937ar','seed',sum(100*clock)));

%% Gather data for all PRNs
% Mega PRN arrays to carry information for all PRNs
OBSTIME_MPRN = zeros(L,32); % In UT
TEC_MPRN = zeros(L,32);
DTEC_MPRN = zeros(L,32);

for kk = 1:1:32,
sv = find(DATA(:,5)==kk);
svdata = DATA(sv,:);

L = length(svdata);
if(isempty(svdata)==0),

    ortw = svdata(:,1); % ort weeks
    orts = svdata(:,2); % ort seconds
    tec = svdata(:,3); % tec  
    dtec = svdata(:,4); % delta tec
    prn = svdata(:,5); %prn
    
    %gps to utc conversion
    [UTC_time, leap_sec, day_of_year]=gps2utc([ortw orts], 0);
    %UTC_time = year month day hour minutes seconds
    UTC_init = UTC_time(1,:);    
    year = UTC_init(1);
    month = UTC_init(2);
    day = UTC_init(3);
    hour = UTC_init(4);
    
    %Observation time seconds (absolute) in UT:
    obstime = (UTC_time(:,4)-hour)*3600 + UTC_time(:,5)*60 + UTC_time(:,6);
        
    NextDay=find(diff(obstime)<0);
    
    if isempty(NextDay)==0,
        obstime = obstime(1:NextDay(1));
        tec = tec(1:NextDay(1));
        dtec = dtec(1:NextDay(1));
    end
    
    L_obs = length(obstime);
    OBSTIME_MPRN(1:L_obs,kk) = obstime;
    TEC_MPRN(1:L_obs,kk) = tec;
    DTEC_MPRN(1:L_obs,kk) = dtec;    
    
    dt = obstime(4)-obstime(3);
    ind_diff= find(diff(obstime)>50);
    disc_ind = [1,ind_diff',length(obstime)];
    Ttime = 0;
    Ttec = 0;
    TPhtec = 0;
    for idf = 1:1:length(disc_ind)-2,
        dTime = obstime(ind_diff(idf)):dt:obstime(ind_diff(idf)+1);
        timei = NaN*dTime';
        teci = ones(size(timei));
        dteci = ones(size(timei));
        
        Ttime = [Ttime;obstime(disc_ind(idf)+1:disc_ind(idf+1));timei];
        Ttec = [Ttec;tec(disc_ind(idf)+1:disc_ind(idf+1));teci];
        
        Totphtec = cumsum(dt*dtec(disc_ind(idf)+1:disc_ind(idf+1)));
        bias = mean(tec(disc_ind(idf)+1:disc_ind(idf+1))-Totphtec);
        
        TPhtec = [TPhtec;Totphtec+bias;dteci]; 
    end

    Ttime = [Ttime;obstime(disc_ind(end-1)+1:disc_ind(end))];
    Ttec = [Ttec;tec(disc_ind(end-1)+1:disc_ind(end))];
    
    Totphtec = cumsum(dt*dtec(disc_ind(end-1)+1:disc_ind(end)));
    bias = mean(tec(disc_ind(end-1)+1:disc_ind(end))-Totphtec);
    TPhtec = [TPhtec;Totphtec+bias];
    
    Ttime = Ttime(Ttime~=0);
    Ttec = Ttec(Ttec~=0);
    TPhtec = TPhtec(TPhtec~=0);

    figure
    plot(Ttime,Ttec,'r');
    hold on
    plot(Ttime,TPhtec,'b');
    %figfontsizes(11,9)
    grid on
    str = strcat('Slant TEC, PRN:',num2str(kk));
    title(str)
    h_legend=legend('Pseudorange TEC',...
    'Phase TEC','Location','NorthEast');
    set(h_legend,'FontSize',fontsz_lg);
    xstring = strcat({'Time [s] after '},num2str(hour),':00 UT on: (mm/dd/yy)',num2str(month),'/'...
        ,num2str(day),'/',num2str(year));
    xlabel(xstring);
    ylabel('TEC [TECU]')

    str_name = strcat('LowRatePHPSTEC_PRN',num2str(kk));
    if(plot_type == 1)
        set(gcf, 'PaperPositionMode', 'auto');
        plotfile = strcat(outdir,str_name,'.eps');
        saveas(gcf,plotfile,'epsc2');
    else
        plotfile = strcat(outdir,str_name,'.png');
        saveas(gca,plotfile);
    end
    close

end

%clear data;
end

toc
% end