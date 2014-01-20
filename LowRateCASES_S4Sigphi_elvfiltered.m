%% LowRateCASES_S4Sigphi.m 
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

% %% Plot low rate plots of S4 and sigma_phi for each PRN. 
% Data are filtered to ignore values below 20 degree elevation angle.

function LowRateCASES_S4Sigphi_elvfiltered(folder_name,signal_type,...
    Receivername)

% %User defined path to the log files of the data.
folder = dir(folder_name);
s={folder.name}; % read all folders and files in the directory

sep = filesep;

%Specify plot type: plot_type = 1 for eps, anything else for png
plot_type = 1;

set(0,'defaultlinelinewidth',1)
% set(0,'defaultaxesfontname','times')
set(0,'defaultaxesfontsize',12)
set(0,'DefaultLineMarkerSize',8);

%%specify signal type
% signal_type = 0; % 0 = L1C, 1 = L2C etc. 
% 0       GPS_L1_CA
% 1       GPS_L2_CM       // GPS L2 civil M code
% 2       GPS_L2_CL
% 3       GPS_L2_CLM
% 4       GPS_L5_I
% 5       GPS_L5_Q
% 6       GPS_L5_IQ       // GPS L5 civil I+Q combined tracking
% 7       GPS_L1_CA_ALT1  // GPS L1 C/A for alternative L1C bank
% 8       CDMA_UHF_PILOT  // Cellular CDMA pilot, I+Q signal
% 9       CDMA_UHF_SYNC   // Cellular CDMA pilot+sync, I+Q signal
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

%Specify Input file names
scintfilename = strcat(folder_name,'scint.log');

txinfofilename = strcat(folder_name,'txinfo.log');

% Only if the file exists, go through next analysis
if exist(scintfilename,'file')
    % Create output folders
    
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
    
    %% Read txinfo data
    txinfo=load(txinfofilename);
    %Column 6 -> health status, 0 is healthy, column 8 -> PRN, column 2 -> 
    %time in seconds 
    
    %% Reading scintillating and reference channels from Scintillation log
    %% file
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
        %This will work well only if there are only 2 repetitions max.
        %ort_sc = ort_sc(ort_sc~=0) will remove the repeated values. But
        %better leave the zeros as they are and recreate the original
        %"scdata" array based on those zeros.
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
        
%         L = length(data);
        if(isempty(data)==0),
            ortw = data(:,1); % ort week
            
            %check for valid files by looking at week number
            %GPS week number cannot be greater than 2000
            %If not valid data, come out
            ind_valid = find(ortw > 2000);
            if(length(ind_valid) > 1),
                string = strcat(s(ii),...
                    ' file has invalid entries for GPS time.');
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
    outfilename = strcat(folder_name,'prn_files_',...
        signal,sep,'Daily_scintdata.mat');
    save(outfilename,'DATA');
    
    clear prn ortw orts lint s4adata s4bdata ...
        s4cdata sphiadata sphibdata sphicdata
    
    for kk = 1:1:32,
        sv = find(DATA(:,10)==kk);
        svdata = DATA(sv,:);
        
        L = length(svdata);
        if(isempty(svdata)==0),
            outfilename = strcat(folder_name,'prn_files_',...
                signal,sep,'scint_prn',num2str(kk),'.mat');
            save(outfilename,'svdata');
            %                 disp(kk);
            
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
            
%             %Observation time seconds after UTC_init hours:
%             obstime = (UTC_time(:,4)-hour)*3600 + ...
%                 UTC_time(:,5)*60 + UTC_time(:,6);
            
%             xstring = strcat({'Time [s] after '},num2str(hour),...
%                 ':00 UT on: (mm/dd/yy)',num2str(month),'/'...
%                 ,num2str(day),'/',num2str(year));
            
            % Time in absolute hours
            obstime = (UTC_time(:,4))+ ...
                UTC_time(:,5)/60 + UTC_time(:,6)/3600;

            xstring = strcat('Time [hr] on: (mm/dd/yy)',num2str(month),'/'...
                ,num2str(day),'/',num2str(year));

            %plotting part of the code
            subplot(2,1,1)
            plot(obstime,s4adata,'r.');
            hold on; plot(obstime,s4bdata,'b+','markersize',4);
            hold on; plot(obstime,s4cdata,'g.');
            grid on
            str = strcat({'Receiver: '},Receivername,{', signal: '},...
                signal,', PRN:',num2str(kk));
            pt=title(str);
            set(pt,'Interpreter', 'none')

            xlabel(xstring);
            ylabel('S_4')
            legend('S_{4_{0-T}}','S_{4_{0-T/2}}','S_{4_{T/2-T}}',...
                'Location','NE')
            
            % Time in absolute hours
            axis([0 24 -inf inf])
            
            subplot(2,1,2)
            plot(obstime,sphiadata,'r.');
            hold on; plot(obstime,sphibdata,'b+','markersize',4);
            hold on; plot(obstime,sphicdata,'g.');
            grid on
            xlabel(xstring);
            ylabel('\sigma_{\phi}')
            legend('\sigma_{\phi_{0-T}}','\sigma_{\phi_{0-T/2}}',...
                '\sigma_{\phi_{T/2-T}}','Location','NE')
            
            % Time in absolute hours
            axis([0 24 -inf inf])
            
            
            %str_name = strcat(cell2mat(s(ii)),'_L1C_PRN',num2str(kk));
%             str_name = strcat('L1C_LowRate1hour_PRN',num2str(kk));
            str_name = strcat('Rx',Receivername,'L1C_LowRate_PRN',num2str(kk));

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
    
else
    disp('scint.log not available')
    
end % if the scint.log file exists, then go through low rate data anlysis.
