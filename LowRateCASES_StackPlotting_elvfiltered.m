%% LowRateCASES_StackPlotting.m 
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

% % Plot low rate stacked plots of S4 and sigma_phi for all PRNs. 
% Data are filtered to ignore values below 20 degree elevation angle.


function LowRateCASES_StackPlotting_elvfiltered(folder_name,signal_type,...
    Receivername)
tic

% %user defined path to the log files of the data.
folder = dir(folder_name);
s={folder.name}; % read all folders and files in the directory

%File separator for the selected operating system
sep = filesep;

%Specify plot type: plot_type = 1 for eps, anything else for png
plot_type = 1;

set(0,'defaultlinelinewidth',1.6)
%set(0,'defaultaxesfontname','times')
% set(0,'defaultaxesfontsize',10)
set(0,'DefaultLineMarkerSize',5);

%%specify signal type
% signal_type = 0; % 0 = L1C, 1 = L2C etc. 
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
fontsz = 12;
marksz = 4;
linestyle= '.';
fontsz_lg = 11; %font size of legends

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
S4ADATA = zeros(1);
S4BDATA = zeros(1);
S4CDATA = zeros(1);
%sigma phi variables
SPHADATA = zeros(1);
SPHBDATA = zeros(1);
SPHCDATA = zeros(1);
LINT = zeros(1);

%Specify Input file names
scintfilename = strcat(folder_name,'scint.log');

txinfofilename = strcat(folder_name,'txinfo.log');

 %% Read txinfo data
    txinfo=load(txinfofilename);
    %Column 6 -> health status, 0 is healthy, column 8 -> PRN, column 2 -> 
    %time in seconds 

%% Reading scintillating and reference channels from Scintillation log file
scfile = load(scintfilename);
%Find data for the particular signal type
scdata = scfile((scfile(:,14)==signal_type),:);
%Corresponding time stamps
ort_sc =(scdata(:,2)+scdata(:,3));

% %Find and remove repeated entries
%1st derivative to find the time slots (wherever there is a huge
%difference, that's the BUG. There are repeated entries inside scin.log
ort_scd = diff(ort_sc);
scdn_ind=find(ort_scd<-20);
for ii = 1:1:length(scdn_ind),
    indii = find(ort_sc == ort_sc(scdn_ind(ii)));
    ort_sc(indii(1):indii(end)-1)=0; 
    %This will work well only if there are only 2 repetitions max. ort_sc =
    %ort_sc(ort_sc~=0) will remove the repeated values. But better leave
    %the zeros as they are and recreate the original "scdata" array based
    %on those zeros.
end
    %Non-zero time stamps means without any repetition.
    nzero_ort = ort_sc~=0;
    scdata = scdata(nzero_ort,:);

% Check for bad values of week number Column # 1 == 0 (must be greater than
% 1600 or so) AND garbage SPR values: Column # 12 (-99, should be greater than
% -80) AND Garbage values = -1 of S4 and sigma_phi (columns 5 to 10)
scdata = scdata((scdata(:,1)>1000 & scdata(:,12)>-80 & scdata(:,5)>-1 ...
    & scdata(:,6)>-1 & scdata(:,7)>-1 & scdata(:,8)>-1 & scdata(:,9)>-1 ...
    & scdata(:,10)>-1),:);


%Find reference satellites
scdatans = scdata((scdata(:,13)==1),:);

ortwa= scdatans(:,1);
ortsa = scdatans(:,2)+scdatans(:,3);
s4a = scdatans(:,5);
sig_phia = scdatans(:,8);

%gps to utc conversion
[UTC_time, leap_sec, day_of_year]=gps2utc([ortwa ortsa], 0);
%UTC_time = year month day hour minutes seconds
UTC_init = UTC_time(1,:);
YEAR = UTC_init(1);
MONTH = UTC_init(2);
DAY = UTC_init(3);
HOUR = UTC_init(4);

obstimea = (UTC_time(:,4)-HOUR)*3600 + UTC_time(:,5)*60 + UTC_time(:,6);
plot(obstimea/3600,s4a,'b.','markersize',10)
xstring = strcat({'Time [hours] after '},num2str(HOUR),...
    ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
        ,num2str(DAY),'/',num2str(YEAR));
xlabel(xstring);
str = strcat({'Low Rate Scintillation Index, S4 of Non-scintillating channels'});
title(str)
axis([0 7 0 0.1])
str_name = strcat(signal,'_LowRateS4_noisefloor_non-scint1');
if(plot_type == 1)
    set(gcf, 'PaperPositionMode', 'auto');
    plotfile = strcat(outdir,str_name,'.eps');
    saveas(gcf,plotfile,'epsc2');
else
    plotfile = strcat(outdir,str_name,'.png');
    saveas(gca,plotfile);
end

figure
plot(obstimea/3600,sig_phia,'g.','markersize',10)
str1 = strcat({'Time [hours] after '},num2str(HOUR),...
    ':00 UT on: (mm/dd/yy)',num2str(MONTH),'/'...
        ,num2str(DAY),'/',num2str(YEAR));
str2 = ' ';
xstring = [str2;str1];
xlabel(xstring);

str = strcat({'Low Rate Scintillation Index, \sigma_{\phi} of Non-scintillating channels'});
title(str)
axis([0 7 0 0.1])
str_name = strcat(signal,'_LowRateSPHI_noisefloor_non-scint1');
if(plot_type == 1)
    set(gcf, 'PaperPositionMode', 'auto');
    plotfile = strcat(outdir,str_name,'.eps');
    saveas(gcf,plotfile,'epsc2');
else
    plotfile = strcat(outdir,str_name,'.png');
    saveas(gca,plotfile);
end


for kk = 1:1:32,
    sv = scdata(:,15)==kk;
    
    %Filter according elevation
        datasv = scdata(sv,:);
        ortsvsec = datasv(:,2);
        
        %filter the values of datasv depending on the elevation(column5) 
        %> 20deg. In txinfo, Column 6 -> health status=0, column 8 -> PRN, 
        % column 2 -> time
        txinfosv=txinfo((txinfo(:,8)==kk & txinfo(:,6)==0),:);
        Lsv = length(datasv);
        
        for iisv = 1:1:Lsv,
            if txinfosv((txinfosv(:,2)==ortsvsec(iisv)),5) < 20,
%                 disp(txinfosv((txinfosv(:,2)==ortsvsec(iisv)),5))
                ortsvsec(iisv)=0;
            end
        end
        
        dataall = datasv(ortsvsec~=0,:);
        
%         %data in one hour
%         data = dataall((dataall(:,2)>filtert1 & ...
%             dataall(:,2)<=filtert1+3600),:);

        %data all
        data = dataall;

%         % make sure that elevation filter worked. 
%         for iisv = 1:1:length(data),
%             if txinfosv((txinfosv(:,2)==data(iisv,2)),5) < 20,
%                 disp(txinfosv((txinfosv(:,2)==data(iisv,2)),5))
%                 data(iisv,2)=0;
%             end
%         end

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

        Lint = data(:,4); % Length of the interval
        s4adata = data(:,5); % S4 values, for the whole interval
        s4bdata = data(:,6); % S4 values, for 1st half of interval
        s4cdata = data(:,7); % S4 values, for 2nd half of interval

        sphiadata = data(:,8); % S4 values, for the whole interval
        sphibdata = data(:,9); % S4 values, for 1st half of interval
        sphicdata = data(:,10); % S4 values, for 2nd half of interval

        prn = data(:,15); %prn

        ORTW = [ORTW, ortw'];
        ORTS = [ORTS, orts'];
        LINT = [LINT, Lint'];
        S4ADATA = [S4ADATA, s4adata'];
        S4BDATA = [S4BDATA, s4bdata'];
        S4CDATA = [S4CDATA, s4cdata'];

        SPHADATA = [SPHADATA, sphiadata'];
        SPHBDATA = [SPHBDATA, sphibdata'];
        SPHCDATA = [SPHCDATA, sphicdata'];

        PRN = [PRN, prn'];

    end
    %
    %                 clear idata qdata rxtime ortw orts prn
    %                 clear indexr rxt_new id_new qd_new prn_new
    %
    %                 %clear data;
end %for loop for prns

%Get rid of the leading zero
ORTW = ORTW(2:end);
ORTS = ORTS(2:end);
LINT = LINT(2:end);
S4ADATA = S4ADATA(2:end);
S4BDATA = S4BDATA(2:end);
S4CDATA = S4CDATA(2:end);
SPHADATA = SPHADATA(2:end);
SPHBDATA = SPHBDATA(2:end);
SPHCDATA = SPHCDATA(2:end);
PRN = PRN(2:end);

DATA = [ORTW; ORTS; LINT; S4ADATA; S4BDATA; S4CDATA; ...
SPHADATA; SPHBDATA; SPHCDATA; PRN]';
% % struct('Data',DATA);
outfilename = strcat(folder_name,'prn_files_',signal,sep,'Daily_scintdata.mat');
save(outfilename,'DATA');

clear prn ortw orts lint s4adata s4bdata s4cdata sphiadata sphibdata sphicdata

ORTS_BEG = ORTS(1);
ORTS_END = ORTS(end);
L = length(DATA);
ORTS_min = min(ORTS);

%gps to utc conversion
[UTC_time, leap_sec, day_of_year]=gps2utc([ORTW(1) ORTS_min], 0);
%UTC_time = year month day hour minutes seconds
UTC_init = UTC_time(1,:);
YEAR = UTC_init(1);
MONTH = UTC_init(2);
DAY = UTC_init(3);
HOUR = UTC_init(4);

%create a proper date string
date = strcat(num2str(MONTH),'/',num2str(DAY),'/',num2str(YEAR));
mmddyy = datestr(datenum(date, 'mm/dd/yy'), 2);

% %to completely randomize the colors of stack plot
% RandStream.setDefaultStream ...
%      (RandStream('mt19937ar','seed',sum(100*clock)));

%% Gather data for all PRNs
% Mega PRN arrays to carry information for all PRNs
OBSTIME_MPRN = zeros(L,32); % In UT
S4A_MPRN = zeros(L,32);
SPHA_MPRN = zeros(L,32);

for kk = 1:1:32,
sv = find(DATA(:,10)==kk);
svdata = DATA(sv,:);

L = length(svdata);
if(isempty(svdata)==0),

    ortw = svdata(:,1); % ort weeks
    orts = svdata(:,2); % ort seconds
    lint = svdata(:,3); % Length of the interval
    
    s4adata = svdata(:,4); % S4 values
    s4bdata = svdata(:,5); % 
    s4cdata = svdata(:,6); % 
    
    sphiadata = svdata(:,7); % Sigma phi values
    sphibdata = svdata(:,8); % 
    sphicdata = svdata(:,9); %
    prn = svdata(:,10); %prn
    
    %gps to utc conversion
    [UTC_time, leap_sec, day_of_year]=gps2utc([ortw orts], 0);
    %UTC_time = year month day hour minutes seconds
    UTC_init = UTC_time(1,:);    
    year = UTC_init(1);
    month = UTC_init(2);
    day = UTC_init(3);
    hour = UTC_init(4);
    
    %Observation time seconds (absolute) in UT:
    obstime = (UTC_time(:,4))*3600 + UTC_time(:,5)*60 + UTC_time(:,6);
    L_obs = length(obstime);
    
    OBSTIME_MPRN(1:L_obs,kk) = obstime;
    S4A_MPRN(1:L_obs,kk) = s4adata;
    SPHA_MPRN(1:L_obs,kk) = sphiadata;    
    
end

%clear data;
end

%find the maximum scintillation that day for nomalization
s4max = max(max(S4A_MPRN));
sphimax = max(max(SPHA_MPRN));

NScale = 4; %the values be scaled by this number

%Displacement vector value on y axis
S4_1 = 1/s4max;
sigmaPhrad1 = 1/sphimax;

%% Plotting Part of the code
count = 1;
figure;
set(gca,'FontSize',fontsz)
for kk = 1:1:32,

    obstime = OBSTIME_MPRN(:,kk);
    s4adata = S4A_MPRN(:,kk);

    obstime = obstime(obstime~=0);
    s4adata = s4adata(s4adata~=0);
    
    [obstime,IX] = sort(obstime);
    s4adata = s4adata(IX);

    if (length(obstime)>=4),
        time = (obstime-HOUR*3600)/3600;
        S4 = (s4adata/s4max)*NScale+count;
        
        dt = time(4)-time(3);
        ind_diff= find(diff(time)>0.5);
        disc_ind = [1,ind_diff',length(time)];
        Ttime = 0;
        SS4 = 0; 
        for idf = 1:1:length(disc_ind)-2,
            dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
            timei = NaN*dTime';
            S4i = ones(size(timei));
            
            Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
            SS4 = [SS4;S4(disc_ind(idf)+1:disc_ind(idf+1));S4i];
        end
        
          Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
          SS4 = [SS4;S4(disc_ind(end-1)+1:disc_ind(end))];  
        
          Ttime = Ttime(Ttime~=0);
          SS4 = SS4(SS4~=0);
          
        plot(Ttime,SS4,...
            'Color',[rand(1) rand(1) rand(1)]);
        PRN_label = (num2str(kk));
        text(time(1), S4(1), PRN_label, 'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right','FontSize',fontsz_lg)
        hold on
        count = count +1;
    end
end

str = strcat({'Low Rate S4 for '},num2str(count-1),...
    {' PRNs'},{'Receiver: '},Receivername,{', signal: '},signal);
            
% str = strcat({'Low Rate Scintillation Index, S4, for '},num2str(count-1),...
%     {' PRNs on '},signal, {' signal'});
pt=title(str);
set(pt,'Interpreter', 'none')

str1 = strcat({'Time [hours] after '},num2str(HOUR),... % [min]
    {':00 UT on (mm/dd/yy) '},mmddyy);
str2 = ' ';
xstring = [str2;str1];
xlabel(xstring);
set(gca,'ytick',[]) 

str1 = {'PRNs'};
str2 = ' ';
ystring = [str1;str2];
ylabel(ystring)
ylabh = get(gca,'YLabel');
set(ylabh,'Position',[-0.1 round(count/2) 0])
%axis([-inf inf 0 count+2])
axis tight

%axis([-inf inf 0 count])
axis tight

%Draw an arrow annotation based on the x and y axis values from the plot
% xPlot = [7.4 7.4];
% yPlot = [15 15+S4_01*NScale];
axPos = get(gca,'Position'); %# gca gets the handle to the current axes

% get the limits, i.e. min and max of the axes.
xMinMax = xlim;
yMinMax = ylim;

%Compute x y values
xPlot = [xMinMax(2)+0.4 xMinMax(2)+0.4];
yPlot = [round(yMinMax(2)/2) round(yMinMax(2)/2)+S4_1*NScale];

% calculate the annotation x and y from the plot x and y.
xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

%txtar = annotation('textarrow',xAn,yAn,'String','S4=1','FontSize',fontsz_lg);


str_name = strcat('Rx',Receivername,signal,'_LR_S4_StackPlot',num2str(count-1),'PRNs');
if(plot_type == 1)
    set(gcf, 'PaperPositionMode', 'auto');
    plotfile = strcat(outdir,str_name,'.eps');
    saveas(gcf,plotfile,'epsc2');
else
    plotfile = strcat(outdir,str_name,'.png');
    saveas(gca,plotfile);
end

%% Plotting Part of the code
count = 1;
figure;
set(gca,'FontSize',fontsz)
for kk = 1:1:32,

    obstime = OBSTIME_MPRN(:,kk);
    sphiadata = SPHA_MPRN(:,kk); 
    
    obstime = obstime(obstime~=0);
    sphiadata = sphiadata(sphiadata~=0);

    [obstime,IX] = sort(obstime);
    sphiadata = sphiadata(IX);
    
    if (length(obstime)>4),
        time = (obstime-HOUR*3600)/3600;
        sigma_PH = sphiadata/sphimax*NScale+count;
        
        dt = time(4)-time(3);
        ind_diff= find(diff(time)>0.5);
        disc_ind = [1,ind_diff',length(time)];
        Ttime = 0;
        Ssigma_PH = 0; 
        for idf = 1:1:length(disc_ind)-2,
            dTime = time(ind_diff(idf)):dt:time(ind_diff(idf)+1);
            timei = NaN*dTime';
            sigma_PHi = ones(size(timei));
            
            Ttime = [Ttime;time(disc_ind(idf)+1:disc_ind(idf+1));timei];
            Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(idf)+1:disc_ind(idf+1));sigma_PHi];
        end
        
          Ttime = [Ttime;time(disc_ind(end-1)+1:disc_ind(end))];
          Ssigma_PH = [Ssigma_PH;sigma_PH(disc_ind(end-1)+1:disc_ind(end))];  
        
          Ttime = Ttime(Ttime~=0);
          Ssigma_PH = Ssigma_PH(Ssigma_PH~=0);
        
        plot(Ttime,Ssigma_PH,...
            'Color',[rand(1) rand(1) rand(1)]);
        PRN_label = (num2str(kk));
        text(time(1), sigma_PH(1), PRN_label, 'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right','FontSize',fontsz_lg)
        hold on
        count = count +1;
    end
end
%Change the fontsize of PRN numbers
    figureHandle = gca;
    set(findall(figureHandle,'type','text'),'fontSize',fontsz_lg)

Rxname = regexp(Receivername, '[_.]', 'split');    
   
Rxnamenew = strcat(Rxname(1),'\_',Rxname(2));

str = strcat({'Low Rate \sigma_{\phi} [normalized] for '},num2str(count-1),...
    {' PRNs'},{'Receiver: '},Rxnamenew,{', signal: '},signal);    
    
% str = strcat({'Low rate standard deviation of phase, \sigma_{\phi} [normalized]'});
%    {' for '},num2str(count-1),...
%    {' PRNs on '},signal,{' signal'});
pt=title(str);
% set(pt,'Interpreter', 'none')

str1 = strcat({'Time [hours] after '},num2str(HOUR),... % [min]
    {':00 UT on (mm/dd/yy) '},mmddyy);
str2 = ' ';
xstring = [str2;str1];
xlabel(xstring);
set(gca,'ytick',[]) 

str1 = {'PRNs'};
str2 = ' ';
ystring = [str1;str2];
ylabel(ystring)
ylabh = get(gca,'YLabel');
set(ylabh,'Position',[-0.2 round(count/2) 0])
%axis([-inf inf 0 count+2])
axis tight

%Draw an arrow annotation based on the x and y axis values from the plot
% xPlot = [7.4 7.4];
% yPlot = [15 15+sigmaPhrad1*NScale];
axPos = get(gca,'Position'); %# gca gets the handle to the current axes

% get the limits, i.e. min and max of the axes.
xMinMax = xlim;
yMinMax = ylim;
%Compute x y values
xPlot = [xMinMax(2)+0.3 xMinMax(2)+0.3];
yPlot = [round(yMinMax(2)/2)-1 round(yMinMax(2)/2)-1+sigmaPhrad1*NScale];

% calculate the annotation x and y from the plot x and y.
xAn = axPos(1) + ((xPlot - xMinMax(1))/(xMinMax(2)-xMinMax(1))) * axPos(3);
yAn = axPos(2) + ((yPlot - yMinMax(1))/(yMinMax(2)-yMinMax(1))) * axPos(4);

%txtar = annotation('line',xAn,yAn);

% %horizontal lines
% hrx1 = [xAn(1)-0.005 xAn(1)+0.005];
% hry1 = [yAn(1) yAn(1)];
% annotation('line',hrx1,hry1);
% 
% hrx2 = [xAn(1)-0.005 xAn(1)+0.005];
% hry2 = [yAn(2) yAn(2)];
% annotation('line',hrx2,hry2);

%text(xMinMax(2)+0.1,round(yMinMax(2)/2)-1.5,'\sigma_{\phi}=1 rad','FontSize',fontsz_lg);
% text(xMinMax(2)+0.1,round(yMinMax(2)/2)-1.1);

str_name = strcat('Rx',Receivername,signal,'_LR_SPHI_StackPlot',num2str(count-1),'PRNs');
if(plot_type == 1)
%     set(gcf, 'PaperPositionMode', 'auto');
    plotfile = strcat(outdir,str_name,'.eps');
%     saveas(gcf,plotfile,'epsc2');
    print('-depsc2','-cmyk',plotfile)
else
    plotfile = strcat(outdir,str_name,'.png');
    saveas(gca,plotfile);
end

toc
