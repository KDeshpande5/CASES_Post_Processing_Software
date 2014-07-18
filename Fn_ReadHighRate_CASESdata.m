%% Fn_ReadHighRate_CASESdata.m 

% This code is a part of CASES_Post_Processing_Software v1
% Copyright (C) 2012 - 2014, Kshitija Deshpande, Virginia Tech. 

% CASES_Post_Processing_Software v1 is free software: you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
% 
% CASES_Post_Processing_Software v1 is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% -----------------------------------------------------------------------

% %% Fn_ReadHighRate_CASESdata.m function reads IQ data first and removes
% repeated entries. For each PRN, considering it is a scintillating
% channel, it finds reference signal for a continuous time, whenever the
% reference channel changes, creates a new segment. Finally, differences
% and saves the data as Filtered .mat files.

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

% % Outputs: 
% PRN_files_<SignalType>/Processed_HR_Data.mat with processed phase, power
% data as a function of time for all PRNs and <SignalType> are saved.


function Fn_ReadHighRate_CASESdata(folder_name,signal_type)
tic

%File separator for the selected operating system
sep = filesep;

%Specify plot type: plot_type = 1 for eps, anything else for png
% plot_type = 0;

%Turn warnings off
warning off  all

%Specify signal type
% signal_type = 0; 
switch signal_type
   case 0
      signal = 'L1CA';
   case 1
      signal = 'L2CM';
   case 2
      signal = 'L2CL';
   case 3
      signal = 'L2CLM';
   case 4
      signal = 'L5I';
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

%Specify path to save plots and other files
% command = strcat('mkdir',{' '},folder_name,'plots_',signal);
% system(cell2mat(command));
% command = strcat('mkdir',{' '},folder_name,'plots_',signal,sep,...
%     'HighRate_ScintAll',sep);
% system(cell2mat(command));

command = strcat('rm -rf',{' '},folder_name,'PRN_files_',signal);
system(cell2mat(command));

command = strcat('mkdir',{' '},folder_name,'PRN_files_',signal);
system(cell2mat(command));
% outdir = strcat(folder_name,'plots_',signal,sep,'HighRate_ScintAll',sep);

%Specify Input file names
iqfilename = strcat(folder_name, 'iq.log');
scintfilename = strcat(folder_name, 'scint.log');

% READ All DATA
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
nzero_ort = find(ort_sc~=0);
scdata = scdata(nzero_ort,:);
%Reread time for scint file(with no repeated entries).
ort_sc = scdata(:,2)+scdata(:,3);

% Check for bad values of week number Column # 1 == 0 (must be greater than
% 1600 or so) AND garbage SPR values: Column # 12 (-99, should be greater than
% -80) AND Garbage values = -1 of S4 and sigma_phi (columns 5 to 10)
scdata = scdata((scdata(:,1)>1000 & scdata(:,12)>-80 & scdata(:,5)>-1 ...
    & scdata(:,6)>-1 & scdata(:,7)>-1 & scdata(:,8)>-1 & scdata(:,9)>-1 ...
    & scdata(:,10)>-1),:);


%% Reading IQ data
tic
iqfile = load(iqfilename);
%Find data for the particular signal type
iqdata = iqfile((iqfile(:,10)== signal_type),:); 
ortw = iqdata(:,3); % ort week
%Check for valid GPS week number entries (cannot be greater than 2000)
ind_valid = find(ortw < 2000);
viqdata = iqdata(ind_valid,:);
orts_iq = viqdata(:,4) + viqdata(:,5); %ort seconds + fractional seconds
toc
% Look inside IQ files for the time slots mentioned in the array of
% structures: ScintStruct for scintillating events, find the scintillating
% satellites during those time and operate on them.

%Initializations
flag_sc = 0;
flag_ref = 0;
count = 0;
REFPRN = zeros(1);
SCINTPRN =zeros(1);
OBSTIME = zeros(1); %UT observation time in seconds
ORTW = zeros(1); %ORT weeks
ORTS = zeros(1); %ORT seconds + fractional seconds
PCDATA = zeros(1); % Processed carrier phase in cycles
PIQPHDATA = zeros(1); % Processed IQ phase data
PIQPOWDATA = zeros(1); % Processed IQ power data

tic
for kk= 1:1:32,
   % Data for scintillating satellite kk
   svdataO = viqdata(viqdata(:,11)==kk,:); %Original (with time repetition)
   
   if ~isempty(svdataO),
   orts = svdataO(:,4) + svdataO(:,5);

   orts_diff = diff(orts);
   orts_ind=find(orts_diff<0);
   for ii = 1:1:length(orts_ind),
       indii = find(orts == orts(orts_ind(ii)));
       orts(indii(1):indii(end)-1)=0;
       %This will work well only if there are only 2 repetitions max. ort_sc =
       %ort_sc(ort_sc~=0) will remove the repeated values. But better leave
       %the zeros as they are and recreate the original "scdata" array based
       %on those zeros.
   end
   nzero_ort = find(orts~=0);
   svdata = svdataO(nzero_ort,:); 
   
   ortw = svdata(:,3);
   ortsSc = svdata(:,4) + svdata(:,5);
   cdataSc = svdata(:,6);
   idata = svdata(:,7);
   qdata = svdata(:,8);
   iqpowSc = idata.^2 + qdata.^2;
   iqphSc = atan2(qdata,idata);

   % Find close-by reference PRNs for kk scintillating PRN from scint.log,
   %also note the respective times of those reference PRNs
   % Low rate times are the end times of the PRNs
   PRNkk_ind = find(scdata(:,13)==0 & scdata(:,15)==kk); 
   %low rate orts
   ortsLR_kk = scdata(PRNkk_ind,2)+ scdata(PRNkk_ind,3);
   % low rate temporary vectors for reference PRNs and times
   reftempLR = find(scdata(:,13)==1);
   ortreftempLR = scdata(reftempLR,2)+scdata(reftempLR,3);
   PRNreftempLR = scdata(reftempLR,15);
   % Find closest reference times to low rate times for kk PRN
%    ortsLR_ref = zeros(size(ortsLR_kk));
%    PRNrefLR = zeros(size(ortsLR_kk));
%    tpind = 0;
%    for jj=1:1:length(ortsLR_kk),
%     [val,indortLR] = min(abs(ortreftempLR-ortsLR_kk(jj)));
%     if val < 20 && val ~=0,
%         tpind = tpind +1;
%         ortsLR_ref(tpind) = ortreftempLR(indortLR);
%         PRNrefLR(tpind) = PRNreftempLR(indortLR);
%     end
%    end
   % % % % M's way
   a = ortreftempLR;
   b = ortsLR_kk;
   [min1,ind1]=min(abs( repmat(a,1,length(b)) - repmat(b,1,length(a))' ));
   ortsLR_ref=a(ind1); %=ortsLR_ref; PRNrefLR = PRNreftempLR(ind1)
   PRNrefLR = PRNreftempLR(ind1);
   
   % Create smaller segments of high rate data depending on the change in
   % reference data
%    a = ortsSc;
%    b = ortsLR_ref;
%    [min1,ind1]=min(abs( repmat(a,1,length(b)) - repmat(b,1,length(a))' ));
%    closeby2b=a(ind1);

   % find continuous data segments in scintillating PRN data
   diff_ind= find(abs(diff(ortsSc))>=1);
   if(~isempty(diff_ind))
       a1 = ortsSc(diff_ind);
       b1 = ortsSc(diff_ind+1);
       endsegkk = [a1; ortsSc(end)];
       stsegkk = [ortsSc(1);b1];
   else
       endsegkk = ortsSc(end);
       stsegkk = ortsSc(1);
   end
   
   % find which reference satellites may be useful for high rate data
   bin_ortsref = zeros(length(stsegkk),length(ortsLR_ref));
   bin_refPRN = zeros(length(stsegkk),length(ortsLR_ref));
   nzero_segkk = zeros(1,length(stsegkk));
   for jj = 1:1:length(stsegkk),
        Ind = find(ortsLR_ref>stsegkk(jj) & ortsLR_ref<endsegkk(jj));
        bin_ortsref(jj,Ind)=ortsLR_ref(Ind);
        bin_refPRN(jj,Ind)=PRNrefLR(Ind);
        if ~isempty(bin_ortsref(jj,bin_ortsref(jj,:)~=0)),
            nzero_segkk(jj) = 1;
%             disp('segment')
%             disp(jj)
%             disp(bin_refPRN(jj,bin_refPRN(jj,:)~=0))
        end
   end
   
   % Segments with some reference channels present based on low rate data
   % Corresponding reference channels in: bin_refPRN and their times:
   % bin_ortsref.
   stsegkkn = stsegkk.*nzero_segkk';
   endsegkkn = endsegkk.*nzero_segkk';
   
   % Inside each segment, look for high rate data corresponding to the
   % reference channel and chop the segment accordingly. Make sure to
   % follow the times from low rate data (e.g. PRN 17 becomes a scint to
   % ref signal from slot 4 to 5 in scint.log file for Jan24)
   indseg = find(stsegkkn~=0);
   
   nseg = length(indseg);
   
   
   for nn = 1:1:nseg,
       
       Rprntp = bin_refPRN(indseg(nn),bin_refPRN(indseg(nn),:)~=0);
       Rorttp = bin_ortsref(indseg(nn),bin_ortsref(indseg(nn),:)~=0);
       dRprntp = diff(Rprntp);
       
       %If there are more than one reference satellites inside on segment,
       %break it into 2
       if(~isempty(dRprntp(dRprntp~=0)))
%            disp('need to break into smaller segments for reference PRNs')
%            disp(Rprntp)
%            disp(kk)

           diffInd_Rprn = find(diff(Rprntp)~=0);
           %number of reference channel segments (may with time discontinuity)
           %switch occurs at Rprntp(diffInd_Rprn)
           for nRefch = 1:1:length(diffInd_Rprn)+1,
               
               if nRefch <= length(diffInd_Rprn),
                    REF = Rprntp(diffInd_Rprn(nRefch));
                    ROrt = Rorttp(diffInd_Rprn(nRefch));
               else
                    REF = Rprntp(diffInd_Rprn(nRefch-1)+1);
                    ROrt = Rorttp(diffInd_Rprn(nRefch-1)+1);
               end
               
               svdataref = viqdata(viqdata(:,11)==REF,:);
               ortsr = svdataref(:,4) + svdataref(:,5);
               
               %find the reference data segment for each PRN
               a = ortsr;               
               if(nRefch == 1)                 
                   b = stsegkkn(indseg(nn));
               else
                   b = Rorttp(diffInd_Rprn(nRefch-1))+0.01;
               end               
               if(nRefch == length(diffInd_Rprn)+1)
                   c = endsegkkn(indseg(nn));
               else
                   c = ROrt;
               end               
               [min1,indst]=min(abs( repmat(a,1,length(b)) - ...
                   repmat(b,1,length(a))' ));
               stsegref=a(indst);
               [min1,indend]=min(abs( repmat(a,1,length(c)) - ...
                   repmat(c,1,length(a))' ));
               endsegref=a(indend);
               
               ortsRef = ortsr(indst:indend);
               cRef = svdataref(indst:indend,6);
               iRef = svdataref(indst:indend,7);
               qRef = svdataref(indst:indend,8);
               iqpowRef = iRef.^2 + qRef.^2;
               iqphRef = atan2(qRef,iRef);
               Lref = length(ortsRef);
               
               if (length(ortsRef)>1) % If no reference data exist, skip
                   
                   % find continuous data segments in scintillating PRN data
                   % corresponding to reference continuous data
                   diff_indr= find(abs(diff(ortsRef))>=1);
                   if(~isempty(diff_indr))
                       % Need to break data (scintillating data according to
                       % discontinuities in reference PRN's time
                       for zmr = 1:1:length(diff_indr)+1,
                           %for all the high rate data
                           %  obstimez = obstime((diff_ind(zm))*(zm-1)+1:...
                           %   zm*diff_ind(zm));
                           if(zmr == 1)                 % first segment
                               zm_st = 1;
                           else
                               zm_st = diff_indr(zmr-1)+1;
                           end
                           
                           if(zmr == length(diff_indr)+1) % last segment
                               zm_end = Lref;
                           else
                               zm_end = diff_indr(zmr);
                           end
                                                      
                           %Assign the segmented data: ref and scint
                           ortsRefz = ortsRef(zm_st:zm_end);
                           cRefz = cRef(zm_st:zm_end);
                           iqpowRefz = iqpowRef(zm_st:zm_end);
                           iqphRefz = iqphRef(zm_st:zm_end);
                           
                           %Get scintillating PRN data close the reference
                           %times
                           a = ortsSc;
                           b = ortsRef(zm_st);
                           c = ortsRef(zm_end);
                           [min1,indst]=min(abs( repmat(a,1,length(b)) ...
                               - repmat(b,1,length(a))' ));
                           %closeb=a(indst);
                           [min1,indend]=min(abs( repmat(a,1,length(c)) ...
                               - repmat(c,1,length(a))' ));
                           %closec=a(indend);
                           ortsScsz = ortsSc(indst:indend);
                           cScsz = cdataSc(indst:indend);
                           iqpowScsz = iqpowSc(indst:indend);
                           iqphScsz = iqphSc(indst:indend);
                           
                           %process only if significant amount of data present
                           if (length(ortsRefz)>10 && length(ortsScsz)>10)
                               
                               %bring the reference and scintillating data to
                               %same size
                               if length(ortsScsz)~=length(ortsRefz),
                                   END = min(length(ortsRefz),...
                                       length(ortsScsz));
                                   cRefz = cRefz(1:END);
                                   cScsz = cScsz(1:END);
                                   iqphRefz = iqphRefz(1:END);
                                   iqphScsz = iqphScsz(1:END);
                                   iqpowRefz = iqpowRefz(1:END);
                                   iqpowScsz = iqpowScsz(1:END);
                               else
                                   END = length(cRefz);
                                   
                               end
                               
                               %Common time
                               Tref_scz = ortsScsz(1:END);
                               
                               %Differenced phases
                               diffcdataz = cScsz - cRefz; 
                               %differenced carrier phase in cycles 
                               
                               diffciqphz = diffcdataz*2*pi + ...
                                   (iqphScsz-iqphRefz);
                               %diff total iq phase in rad
                               
                               %Raw power
                               rawiqpowz = iqpowScsz;
                               
                               scint_prnz = kk * ones(size(Tref_scz));
                               ref_prnz = Rprntp(1) * ones(size(Tref_scz));
                               ortwscz = ortw(1) * ones(size(Tref_scz));
 
                               %Append the processed data to arrays to
                               %later plot for each PRN
                               ORTW = [ORTW, ortwscz'];
                               ORTS = [ORTS, Tref_scz'];
                               PCDATA = [PCDATA, diffcdataz'];
                               PIQPHDATA = [PIQPHDATA, diffciqphz'];
                               PIQPOWDATA = [PIQPOWDATA, ...
                                   rawiqpowz'];
                               SCINTPRN = [SCINTPRN, scint_prnz'];
                               REFPRN = [REFPRN, ref_prnz'];
                           end % if length of segment is smaller than 10 samples
                       end % for all continuous time segments inside reference data
                       
                   else % if no discontinuities in data for this ref PRN
                       
                       %Get scintillating PRN data close the reference times
                       a = ortsSc;
                       b = ortsRef(1);
                       c = ortsRef(end);
%                        a = ortsSc;
%                        b = stsegkkn(indseg(nn));
%                        c = endsegkkn(indseg(nn));
                       [min1,indst]=min(abs( repmat(a,1,length(b)) ...
                           - repmat(b,1,length(a))' ));
                       closeb=a(indst);
                       [min1,indend]=min(abs( repmat(a,1,length(c)) ...
                           - repmat(c,1,length(a))' ));
                       closec=a(indend);
                       ortsScs = ortsSc(indst:indend);
                       cScs = cdataSc(indst:indend);
                       iqpowScs = iqpowSc(indst:indend);
                       iqphScs = iqphSc(indst:indend);
                       
                       %bring the reference and scintillating data to same size
                       if length(ortsScs)~=length(ortsRef),
                           END = min(length(ortsRef),length(ortsScs));
                           cRef = cRef(1:END);
                           cScs = cScs(1:END);
                           iqphRef = iqphRef(1:END);
                           iqphScs = iqphScs(1:END);
                           iqpowRef = iqpowRef(1:END);
                           iqpowScs = iqpowScs(1:END);
                       else
                           END = length(cRef);
                           
                       end
                       
                       %Common time
                       Tref_sc = ortsScs(1:END);
                       
                       %Differenced phases
                       diffcdata = cScs - cRef;
                       
                       diffciqph = diffcdata*2*pi + (iqphScs-iqphRef);
                       
                       %Raw power
                       rawiqpow = iqpowScs;
                       
                       scint_prn = kk * ones(size(Tref_sc));
                       ref_prn = Rprntp(1) * ones(size(Tref_sc));
                       ortwsc = ortw(1) * ones(size(Tref_sc));
                       %                 figure
                       %                 plot(Tref_sc,diffcdata,'r.')
                       %                 title(num2str(kk))
                       % %                 hold on
                       % %                 pause
                       % %                 plot(Tref_sc,cRef,'r.')
                       %                 pause
                       
                       %Append the processed data to arrays to
                       %later plot for each PRN
                       ORTW = [ORTW, ortwsc'];
                       ORTS = [ORTS, Tref_sc'];
                       PCDATA = [PCDATA, diffcdata'];
                       PIQPHDATA = [PIQPHDATA, diffciqph'];
                       PIQPOWDATA = [PIQPOWDATA, ...
                           rawiqpow'];
                       SCINTPRN = [SCINTPRN, scint_prn'];
                       REFPRN = [REFPRN, ref_prn'];
                   end % If smaller segments in reference data (due to
                   % discontinuities in reference times)
               end % If no reference data exist, skip
               
               
               
           end % for: nRefch: number of different reference signals inside
           
           
           
       else % Just one reference PRN
            svdataref = viqdata(viqdata(:,11)==Rprntp(1),:);
            ortsr = svdataref(:,4) + svdataref(:,5);
            
            a = ortsr;
            b = stsegkkn(indseg(nn));
            c = endsegkkn(indseg(nn));
            [min1,indst]=min(abs( repmat(a,1,length(b)) - ...
                repmat(b,1,length(a))' ));
            stsegref=a(indst);
            [min1,indend]=min(abs( repmat(a,1,length(c)) - ... 
                repmat(c,1,length(a))' ));
            endsegref=a(indend);
            
            ortsRef = ortsr(indst:indend);
            cRef = svdataref(indst:indend,6);
            iRef = svdataref(indst:indend,7);
            qRef = svdataref(indst:indend,8);
            iqpowRef = iRef.^2 + qRef.^2;
            iqphRef = atan2(qRef,iRef);
            Lref = length(ortsRef);
            
            if (length(ortsRef)>10) % If no reference data exist, skip
                
                % find continuous data segments in scintillating PRN data
                % corresponding to reference continuous data
                diff_indr= find(abs(diff(ortsRef))>=1);
                if(~isempty(diff_indr))
                   % Need to break data (scintillating data according to
                   % discontinuities in reference PRN's time
                   for zmr = 1:1:length(diff_indr)+1,
                       %for all the high rate data
                       %  obstimez = obstime((diff_ind(zm))*(zm-1)+1:...
                       %   zm*diff_ind(zm));
                       if(zmr == 1)                 % first segment
                           zm_st = 1;
                       else
                           zm_st = diff_indr(zmr-1)+1;
                       end
                       
                       if(zmr == length(diff_indr)+1) % last segment
                           zm_end = Lref;
                       else
                           zm_end = diff_indr(zmr);
                       end
                       
                       
                       %Assign the segmented data: ref and scint
                       ortsRefz = ortsRef(zm_st:zm_end);
                       cRefz = cRef(zm_st:zm_end);
                       iqpowRefz = iqpowRef(zm_st:zm_end);
                       iqphRefz = iqphRef(zm_st:zm_end);                       
                       
                       %Get scintillating PRN data close the reference
                       %times
                       a = ortsSc;
                       b = ortsRef(zm_st);
                       c = ortsRef(zm_end);
                       [min1,indst]=min(abs( repmat(a,1,length(b)) ...
                           - repmat(b,1,length(a))' ));
                       %closeb=a(indst);
                       [min1,indend]=min(abs( repmat(a,1,length(c)) ...
                           - repmat(c,1,length(a))' ));
                       %closec=a(indend);
                       ortsScsz = ortsSc(indst:indend);
                       cScsz = cdataSc(indst:indend);
                       iqpowScsz = iqpowSc(indst:indend);
                       iqphScsz = iqphSc(indst:indend);
                       
                       %process only if significant amount of data present
                       if (length(ortsRefz)>10 && length(ortsScsz)>10)
                           
                           %bring the reference and scintillating data to
                           %same size
                           if length(ortsScsz)~=length(ortsRefz),
                               END = min(length(ortsRefz),...
                                   length(ortsScsz));
                               cRefz = cRefz(1:END);
                               cScsz = cScsz(1:END);
                               iqphRefz = iqphRefz(1:END);
                               iqphScsz = iqphScsz(1:END);
                               iqpowRefz = iqpowRefz(1:END);
                               iqpowScsz = iqpowScsz(1:END);
                           else
                               END = length(cRefz);
                               
                           end
                           
                           %Common time
                           Tref_scz = ortsScsz(1:END);
                           
                           %Differenced phases
                           diffcdataz = cScsz - cRefz;
                           
                           diffciqphz = diffcdataz*2*pi + ...
                               (iqphScsz-iqphRefz);
                           
                           %Raw power
                           rawiqpowz = iqpowScsz;
                           
                           scint_prnz = kk * ones(size(Tref_scz));
                           ref_prnz = Rprntp(1) * ones(size(Tref_scz));
                           ortwscz = ortw(1) * ones(size(Tref_scz));
                           %                 figure
                           %                 plot(Tref_sc,diffcdata,'r.')
                           %                 title(num2str(kk))
                           % %                 hold on
                           % %                 pause
                           % %                 plot(Tref_sc,cRef,'r.')
                           %                 pause
                           
                           %Append the processed data to arrays to
                           %later plot for each PRN
                           ORTW = [ORTW, ortwscz'];
                           ORTS = [ORTS, Tref_scz'];
                           PCDATA = [PCDATA, diffcdataz'];
                           PIQPHDATA = [PIQPHDATA, diffciqphz'];
                           PIQPOWDATA = [PIQPOWDATA, ...
                               rawiqpowz'];
                           SCINTPRN = [SCINTPRN, scint_prnz'];
                           REFPRN = [REFPRN, ref_prnz'];
                       end % if length of segment is smaller than 10 samples
                   end % for all continuous time segments inside reference data 
                    
                else
                    % If no discontinuities in the reference PRN data
                    %Get scintillating PRN data close the reference times
                    a = ortsSc;
                    b = stsegkkn(indseg(nn));
                    c = endsegkkn(indseg(nn));
                    [min1,indst]=min(abs( repmat(a,1,length(b)) ...
                        - repmat(b,1,length(a))' ));
                    closeb=a(indst);
                    [min1,indend]=min(abs( repmat(a,1,length(c)) ...
                        - repmat(c,1,length(a))' ));
                    closec=a(indend);
                    ortsScs = ortsSc(indst:indend);
                    cScs = cdataSc(indst:indend);
                    iqpowScs = iqpowSc(indst:indend);
                    iqphScs = iqphSc(indst:indend);
                    
                    %bring the reference and scintillating data to same size
                    if length(ortsScs)~=length(ortsRef),
                        END = min(length(ortsRef),length(ortsScs));
                        cRef = cRef(1:END);
                        cScs = cScs(1:END);
                        iqphRef = iqphRef(1:END);
                        iqphScs = iqphScs(1:END);
                        iqpowRef = iqpowRef(1:END);
                        iqpowScs = iqpowScs(1:END);
                    else
                        END = length(cRef);
                        
                    end
                    
                    %Common time
                    Tref_sc = ortsScs(1:END);
                    
                    %Differenced phases
                    diffcdata = cScs - cRef;
                    
                    diffciqph = diffcdata*2*pi + (iqphScs-iqphRef);
                    
                    %Raw power
                    rawiqpow = iqpowScs;
                    
                    scint_prn = kk * ones(size(Tref_sc));
                    ref_prn = Rprntp(1) * ones(size(Tref_sc));
                    ortwsc = ortw(1) * ones(size(Tref_sc));

                    %Append the processed data to arrays to
                    %later plot for each PRN
                    ORTW = [ORTW, ortwsc'];
                    ORTS = [ORTS, Tref_sc'];
                    PCDATA = [PCDATA, diffcdata'];
                    PIQPHDATA = [PIQPHDATA, diffciqph'];
                    PIQPOWDATA = [PIQPOWDATA, ...
                        rawiqpow'];
                    SCINTPRN = [SCINTPRN, scint_prn'];
                    REFPRN = [REFPRN, ref_prn'];
                end % If smaller segments in reference data (due to
                % discontinuities in reference times)
            end % If no reference data exist, skip
       
       end % for smaller segments with different reference PRN

   end% for all bigger segments       
   
   end% only if there exists any high rate data for satellite kk
end

%Get rid of the leading zero
ORTW = ORTW(2:end);
ORTS = ORTS(2:end);
PCDATA = PCDATA(2:end);
PIQPHDATA = PIQPHDATA(2:end);
PIQPOWDATA = PIQPOWDATA(2:end);
SCINTPRN = SCINTPRN(2:end);
REFPRN = REFPRN(2:end);
OBSTIME = ones(size(ORTS));

DATA = [OBSTIME; ORTW; ORTS; PCDATA; PIQPHDATA; PIQPOWDATA; REFPRN;...
    SCINTPRN]';
outfilename = strcat(folder_name,'PRN_files_',signal,sep,...
    'Processed_HR_Data.mat');
save(outfilename,'DATA');


toc