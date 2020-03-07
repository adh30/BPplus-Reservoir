% BPplus_v12
%% batch analysis of BP+ to give results like Sphygmocor
%  and do reservoir analysis using kreservoir_vXX 
%% Copyright 2019 Alun Hughes based on original code by Kim Parker
% Also uses xml2struct.m by W. Falkena, ASTI, TUDelft, 21-08-2010 with additional 
% modifications by A. Wanner, I Smirnov & Chao-Yuan Yeh
%
% This software is distributed under under the terms of the GNU General Public License
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% http://www.gnu.org/licenses/gpl.html

%% Versions
%  v0.1 (03/01/20) First solid version but still in early beta - no WIA, no
%  save
%  v0.11 (03/01/20) version with WIA but no save except figures
%  v0.12 (03/01/20) save results in excel file and achieve concordance
%  with batch reservoir v13 (sphygmocor) except for HRV
%  v0.12 (07/03/20) add SEVR [Buckberg index] calculation and add quality
%  index to output
%%%%%%%%%%%%%%%% 
%% m files required to be in directory
% ai_v1
% xml2struct
% fitres_v5
% kreservoir_v13
%%%%%%%%%%%%%%%%
%% Constants
    kres_v='v13';          % Version tracking for reservoir fitting
    headernumber=53;       % headers for columns of results (see end)
    mmHgPa = 133;          % P conversion for WIA
    uconst=1;              % empirical constant to convert normalized velocity to m/s
    
%% Select files
% folder_name = uigetdir; % dropped for now as it slows the process
folder_name ='C:\BPPdata\'; % speeds things up
file_lists=dir(fullfile(folder_name, '*.xml'));
no_of_files=length(file_lists);
% set record number to 1 and extract filename
record_no=1;
% preallocate cell array
proc_var=cell(no_of_files,headernumber);

%% Start loop through all files
for file_number=1:no_of_files
% refresh filename 
filename=file_lists(record_no).name;
[data] = xml2struct([folder_name filename]);

% progress bar
    % open waitbar
    if record_no==1
        h = waitbar(record_no/no_of_files,'Processing files...');
    else
        waitbar(record_no/no_of_files,h,'Processing files...');
    end

% extract values
ba_sbp=str2double(data.CardioScope.MeasDataLogger.Sys.Text);                % brachial systolic BP, mmHg
dbp=str2double(data.CardioScope.MeasDataLogger.Dia.Text);                   % diastolic BP, mmHg
ba_pp=ba_sbp-dbp;                                                           % brachial pulse pressure, mmHg
map=str2double(data.CardioScope.MeasDataLogger.Mean.Text);                  % mean arterial pressure, mmHg 
hr=str2double(data.CardioScope.MeasDataLogger.Hr.Text);                     % heart rate, bpm
samplerate=str2double(data.CardioScope.MeasDataLogger.SampleRate.Text);     % sample rate, Hz
aosbp=str2double(data.CardioScope.Results.Result.aoSys.Text);               % cSBP calculated by BP+, mmHg
aopp=aosbp-dbp;                                                             % cPP calculated by BP+, mmHg
snr=str2double(data.CardioScope.Results.Result.SNR.Text);                   % Signal to noise ratio, dB
ss_rmssd=str2double(data.CardioScope.Results.Result.RMSSD.Text);            % RMSSD from suprasystolic signal
ss_ai=str2double(data.CardioScope.Results.Result.ssAI.Text);                % AI from suprasystolic signal
ssdpdt=str2double(data.CardioScope.Results.Result.ssDpDtMax.Text);          % dp/dt from suprasystolic signal in uncorrected units
ssHARM=str2double(data.CardioScope.Results.Result.ssHARM.Text);             % unknown 
ssPP=str2double(data.CardioScope.Results.Result.ssPP.Text);                 % unknown
ssPPV=str2double(data.CardioScope.Results.Result.ssPPV.Text);               % Pulse pressure variation, % [?]
ssRWTTFoot=str2double(data.CardioScope.Results.Result.ssRWTTFoot.Text);     % reflected wave transit time from foot of suprasystolic signal
ssRWTTPeak=str2double(data.CardioScope.Results.Result.ssRWTTPeak.Text);     % reflected wave transit time from peak of suprasystolic signal
ssSEP=str2double(data.CardioScope.Results.Result.ssSEP.Text);               % Systolic ejection period, s
algo=data.CardioScope.Results.Result.Attributes.algorithm_revision;         % Software algorithm
ssTn=split(data.CardioScope.Results.Result.ssAverageBeatPointsIdxs.Text,','); % Times of characteristic points?
% Timings of peaks dont match waveforms - seems to be a bug or possibly
% timings are from foot rather than beginning of waveform?
T1=str2double(ssTn(2));                                                     % Time of 1st peak in samples
T2=str2double(ssTn(3));                                                     % Time of inflection in samples
T3=str2double(ssTn(4));                                                     % Time of 2nd peak in samples
T4=str2double(ssTn(5));                                                     % Time of nadir of dichrotic notch
datestring=data.CardioScope.MeasDataLogger.Attributes.datetime;             % Date as text string
%% categorise quality based on SNR
if snr>=12
    quality='Excellent';
elseif snr>=9
    quality='Good';
elseif snr>=6
    quality='Acceptable';
elseif snr>0
    quality='Poor';
else 
    elsequality='Unacceptable';
end

%% brachial average beat
sstxt = split(data.CardioScope.Results.Result.ssAverageBeat.Text,',');
%b_avp_av=zeros(size(sstxt,1),size(sstxt,2)); 
ba_p_av=str2double(sstxt);
calss_p=ba_pp/(max(ba_p_av)-min(ba_p_av));
ba_p_av=dbp+(ba_p_av*calss_p);

ssdpdt=ssdpdt*calss_p;                                                      % correcting to mmHg/s

%% aortic average beat, ao_p_av
sstxt = split(data.CardioScope.Results.Result.aoAverageBeat.Text,',');
%b_avp_av=zeros(size(sstxt,1),size(sstxt,2)); 
a=str2double(sstxt);

% fix disjoin
[~, newend]=min(a);
zx=length(a)-newend;
extra=3;
aa=nan(length(a)+extra,1);      % add an extra 3 nan elements to allow for the crop
aa(1:zx)=a(newend+1:end);
aa(zx+extra+1:length(aa))=a(1:newend);
ao_p_av = fillmissing(aa,'spline')';      % use fillmissing to fill with a spline

%% calculate aoAIx
[~, yy]= max((gradient(gradient(ao_p_av))));    % define start as peak 
%dp^2 since some of the aortic pressures have a hump at the start
% p is variable to pass
p=ao_p_av(yy:end);
[ao_ai, aoPi, aoTfoot, aoTi, aoTmax, aoTypetxt] = ai_v1(p, samplerate);

%% calculate SEVR (aortic and brachial)
ao_sevr=trapz(ao_p_av(lsys+1:end))/trapz(ao_p_av(1:lsys));
ba_sevr=trapz(ba_p_av(lsys+1:end))/trapz(ba_p_av(1:lsys));

%% calculate some more measures
[~, Tpeaks]=findpeaks(diff(p), 'Npeaks', 3);    % peaks of dP
aosbp2=p(Tpeaks(2));          % this is more reliable than 3rd zero crossing of 4th derivative
%[~,aoTd]=min(diff(p));       % MAY NOT BE NECESSARY TO CALCULATE min dP/dt TWICE 
% ao_es=p(aoTd);
ao_dpdt=max(diff(p))*samplerate;

%% do reservoir calculations for aortic pressure
[aoTn_av, aoPinf_av, aoP_av,aoPr_av,aoPn_av, aofita_av,...
aofitb_av, aorsq_av]=BPfitres_v1(p,samplerate);  
aoPxs=aoP_av-aoPr_av;

%% do AIx reservoir calculations for brachial pressure
% assign p to pass to function
p=ba_p_av';
[baTn_av, baPinf_av, baP_av,baPr_av,baPn_av, bafita_av,...
bafitb_av, barsq_av]=BPfitres_v1(p,samplerate);  
baPxs=baP_av-baPr_av;

% pAI
[~, Tpeaks]=findpeaks(diff(p), 'Npeaks', 3);    % peaks of dP
ba_sbp2=p(Tpeaks(2));
ba_ai=100*(ba_sbp2-dbp)/(ba_sbp-dbp); % peripheral augmentation index, %
%[~,ba_Td]=min(diff(p));            % Duplicate
% ba_es=p(ba_Td);
ba_dpdt=max(diff(p))*samplerate;

%% Create estimates of Pf and Pb using the assumption that Pb = cPr/2
aoPb_av=(aoPr_av-min(aoP_av))/2;
aoPf_av=aoP_av-aoPb_av-min(aoP_av);
Pb_Pf=max(aoPb_av)/max(aoPf_av);
RI=max(aoPb_av)/(max(aoPf_av)+ max(aoPb_av));

%% wave intensity analysis (only done on central P)
% currently just using gradients - may be worth using SG
% correct units 1mmHg = 133Pa to calculate W/m2
% a bit of tweaking to get WIA to work given the windows for the
% derivatives - use a 9 sample shift should work for fsg721
p=aoP_av;
p(9:end)=aoP_av(1:end-8);
p(1)=min(aoP_av);
p(2:8)= nan;
p=fillmissing(p,'spline')'; 

% dp=gradient(p*mmHgPa);
dp=fsg721(p');
% convert Pxs to flow velocity and assume peak velocity = 1m/s
% based on Lindroos, M., et al. J Am Coll Cardiol 1993; 21: 1220-1225.
cu=uconst*(aoPxs/max(aoPxs));
u=cu;
u(9:end)=cu(1:end-8);
u(1:3)=0;
u(4:9)= nan;
u=fillmissing(u,'spline')';

% duxs=gradient(cu);
duxs=fsg721(u');
% check alignment duxs, dp
[~, x1]=max(dp);
[~, x2]=max(duxs);
sh=x1-x2;
duxs(sh+1:end)=duxs(1:end-sh);
di=dp.*duxs;
minpeak=max(di)/10;
    % I've left the warning when there are no peaks for now but if it
    % needs to be suppressed then 'signal:findpeaks:largeMinPeakHeight' is
    % its id
    [dippks,diplocs,dipw]=findpeaks(di, 'NPeaks',2,'MinPeakHeight',minpeak); % find two dIpositive peaks
    [dimpks,dimlocs,dimw]=findpeaks(-di, 'NPeaks',1,'MinPeakHeight',0.7*max(-di)); % find one dIpositive peaks
    di=di*mmHgPa*length(dp);% units fixed - now in W/m2 per cycle^2
    dippks=dippks*mmHgPa*length(dp);
    dimpks=dimpks*mmHgPa*length(dp);
    dipt=diplocs/samplerate;
    dimt=dimlocs/samplerate;
    % calculate areas
    % For a Gaussian curve (assumed) the area is 1.06447*height*width
    diparea=1.06447*dippks.*dipw;
    dimarea=1.06447*dimpks.*dimw;
    wri=dimarea/diparea(1);
    % error trap when W2 is unmeasureable
     if length(dippks)==1
        dippks(2)=0;
        dipt(2)=0;
        diparea(2)=0;
     end
    
    
    % Estimate c (wavespeed) as k*dP/du where k is empirical constant
    % currently k = 1!
    rhoc=max(aoPxs)*mmHgPa/1000; % fixed units (m/s)
%% possible problems with fits
    if aoPinf_av>=dbp || aoPinf_av<-12 || aofita_av<=0 || aofitb_av<=0 || aorsq_av<0.9
        prob=1;
    else
        prob =0;
    end      
% *************should write error salvage sometime!

%% make figures and data subfolders
    figfolder='C:\BPPdata\figures\';
    datafolder='C:\BPPdata\results\'; 
    if ~exist(figfolder, 'dir')
    mkdir(figfolder);
    end
    if ~exist(datafolder, 'dir')
    mkdir(datafolder);
    end
    % save separate jpg (for viewing) and wmf (for editing)
    expression = 'xml';
    replace1 = 'wmf';
    replace2 ='jpg';
    wmffile = regexprep(filename,expression,replace1);  
    jpgfile = regexprep(filename,expression,replace2);   
    
%% Print figures
    % make a time variable for printing
    Time=(1:length(aoP_av))/samplerate;
    % Aortic pressure and reservoir
    f1=figure('visible','off');                     % dont display figure
    plot(Time,aoP_av-aoP_av(1),Time,aoPr_av-aoPr_av(1),'r', Time, aoPxs,'k')  
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    box off;
    print ('-dmeta', '-r300' , [figfolder wmffile]);
    %print ('-djpeg', '-r300' , [figfolder jpgfile]);
    
    % WI
    f2=figure('visible','off');                     % dont display figure
    TimeDI=Time(1:length(di));
    plot (TimeDI, di);                              % **** to allow for new length
    xlabel('Time (s)')
    ylabel('dI (W/m2)')
    box off;
    wmffile1 = regexprep(filename,'.xml','w.wmf');
    %jpgfile1 = regexprep(filename,'.xml','w.jpg');
    print ('-dmeta', '-r300' , [figfolder wmffile1]);
    % print ('-djpeg', '-r300' , [figfolder jpgfile1]);
    drawnow();                  % added to attempt to stop java leak
    
    % P, Pf, Pb
    f3=figure('visible','off');                     % dont display figure
    hold on; plot(Time, aoPf_av,'b', Time, aoPb_av,'r', Time, aoP_av-min(aoP_av),'k');
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    wmffile2 = regexprep(filename,'.xml','fb.wmf');
    print ('-dmeta', '-r300' , [figfolder wmffile2]);
    drawnow();                  % added to attempt to stop java leak
    
    % clear and close figures
    clear f1 f2
%     clear f1 f2 f3
    figs =  findobj('type','figure');
    close(figs);
    clear figs;
%% Save results
[~,max_tpb]= max(baP_av);
[~,max_tpa]= max(aoP_av);
[maxPra,max_tra]= max(aoPr_av);
[maxPxsa,max_txsa]= max(aoPxs);
[maxPrb,max_trb]= max(baPr_av);
[maxPxsb,max_txsb]= max(baPxs);

% write variables 
    proc_var{record_no,1}=filename;                                         % filename
    proc_var{record_no,2}=ba_sbp;                                           % brachial SBP, mmHg
    proc_var{record_no,3}=max_tpb/samplerate;                               % time max brachial P (baSBP), s
    proc_var{record_no,4}=dbp;                                              % min P, DBP, mmHg
    proc_var{record_no,5}=sum(aoPr_av)/samplerate;                          % integral aoPres, mmHg.s
    proc_var{record_no,6}=maxPra;                                           % max aoPres (aSBP), mmHg
    proc_var{record_no,7}=max_tra/samplerate;                               % Time max aoPres, sec
    proc_var{record_no,8}=sum(aoPr_av-dbp)/samplerate;                      % integral aoPres-diastolic, mmHg.s
    proc_var{record_no,9}=datestring;                                       % Date as text string
    proc_var{record_no,10}=samplerate;                                      % sampling rate, 1/sec [used to be time of max Pres-diastolic but since this == Time max Pres I replaced it with the sampling rate]
    proc_var{record_no,11}=sum(aoPxs)/samplerate;                           % Integral excess pressure, mmHg.s
    proc_var{record_no,12}=maxPxsa;                                         % max excess P, mmHg
    proc_var{record_no,13}=max_txsa/samplerate;                             % time of max excess P, sec
    proc_var{record_no,14}=aoTn_av;                                         % Time of end systole by max -dp/dt, sec
    proc_var{record_no,15}=aoPinf_av;                                       % aortic P infinity, mmHg
    proc_var{record_no,16}=aoPn_av;                                         % aortic P at end systole by max -dp/dt, mmHg
    proc_var{record_no,17}=aofita_av;                                       % aortic rate constant A (ka), 1/sec
    proc_var{record_no,18}=aofitb_av;                                       % aortic rate constant B (kb), 1/sec
    proc_var{record_no,19}=aorsq_av;                                        % aortic R2 for diastolic fit
    proc_var{record_no,20}=prob;                                            % likely problem with fit
    proc_var{record_no,21}=kres_v;                                          % version tracking
    proc_var{record_no,22}=aoTypetxt;                                       % AI type
    proc_var{record_no,23}=hr;                                              % Heart rate, bpm
    proc_var{record_no,24}=ba_sbp2;                                         % SBP2, mmHg
    proc_var{record_no,25}=sum(baPr_av)/samplerate;                         % integral baPres, mmHg.s
    proc_var{record_no,26}=maxPrb;                                          % max baPres, mmHg
    proc_var{record_no,27}=max_trb/samplerate;                              % Time max Pres, sec
    proc_var{record_no,28}=sum(baPxs)/samplerate;                           % Integral ba excess pressure, mmHg.s
    proc_var{record_no,29}=maxPxsb;                                         % max ba excess P, mmHg
    proc_var{record_no,30}=max_txsb/samplerate;                             % time of max ba excess P, sec
    proc_var{record_no,31}=bafita_av;                                       % brachial rate constant A (ka), 1/sec
    proc_var{record_no,32}=bafitb_av;                                       % brachial rate constant B (kb), 1/sec
    proc_var{record_no,33}=barsq_av;                                        % brachial R2 for diastolic fit
    proc_var{record_no,34}=baPinf_av;                                       % brachial Pinf, mmHg
    proc_var{record_no,35}=baPn_av;                                         % brachial P at end systole by max -dp/dt, mmHg
    proc_var{record_no,36}=ao_dpdt;                                         % aortic dp/dt, mmHg/s
    proc_var{record_no,37}=ba_dpdt;                                         % brachial dp/dt mmHg/s
    proc_var{record_no,38}=Pb_Pf;                                           % Pb/Pf (ratio of backward to forward pressure - also know as reflection magnitude. RM
    proc_var{record_no,39}=RI;                                              % Reflection index (RI) calculated as Pb/(Pb+Pf)
    proc_var{record_no,40}=dippks(1);                                       % W1 intensity
    proc_var{record_no,41}=dipt(1);                                         % W1 time
    proc_var{record_no,42}=diparea(1);                                      % W1 area
    proc_var{record_no,43}=dimpks;                                          % W-1 intensity
    proc_var{record_no,44}=dimt;                                            % W-1 time
    proc_var{record_no,45}=dimarea;                                         % W-1 area
    proc_var{record_no,46}=dippks(2);                                       % W2 intensity
    proc_var{record_no,47}=dipt(2);                                         % W2 time
    proc_var{record_no,48}=diparea(2);                                      % W2 area
    proc_var{record_no,49}=wri;                                             % WRI
    proc_var{record_no,50}=rhoc;                                            % rhoc
    proc_var{record_no,51}=ao_sevr;                                         % aortic SEVR
    proc_var{record_no,52}=ba_sevr;                                         % brachial SEVR
    proc_var{record_no,53}=quality;                                         % quality index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% increment record number
    if record_no==no_of_files
       close(h)
       clear h; 
    else
       record_no=record_no+1;
    end
%% Save the results as an excel spreadheet
% % % % 
% 1 're_file' 
% 2 're_basbp' 
% 3 're_tbasbp' 
% 4 're_minp' 
% 5 're_intaopr' 
% 6 're_maxaopr'
% 7 're_tmaxaopr' 
% 8 're_intaoprlessdbp' 
% 9 'date' 
% 10 're_sam_rate'
% 11 're_intaoxsp' 
% 12 're_maxaoxsp' 
% 13 're_tmaxaoxsp'
% 14 're_aotn' 
% 15 're_aopinf' 
% 16 're_aopn'
% 17 're_aofita'
% 18 're_aofitb'
% 19 're_aorsq'
% 20 're_prob'
% 21 're_version'
% 22 're_aitype' 
% 23 're_hr' 
% 24 're_sbp2' 
% 25 're_intbapr'
% 26 're_maxbapr'
% 27 're_tmaxbapr'
% 28 're_intbaxsp'
% 29 're_maxbaxsp'
% 30 're_tmaxbap'
% 31 're_bafita'
% 32 're_bafitb'
% 33 're_barsq'
% 34 're_bapinf' 
% 35 're_bapn' 
% 36 're_ao_dpdt'
% 37 're_ba_dpdt'
% 38 're_pb_pf' 
% 39 're_ri' 
% 40 're_wf1i' 
% 41 're_wf1t' 
% 42 're_wf1a'
% 43 're_wbi' 
% 44 're_wbt' 
% 45 're_wba' 
% 46 're_wf2i' 
% 47 're_wf2t' 
% 48 're_wf2a'
% 49 're_wri'
% 50 're_rhoc'
% 51 'ao_sevr'
% 52 'ba_sevr'
% 53 'quality'
% % % % 
xlsfile='C:\BPPdata\results\resdata.xls';
header = {'re_file' 're_basbp' 're_tbasbp' 're_minp' 're_intaopr' 're_maxaopr'...
    're_tmaxaopr' 're_intaoprlessdbp' 'date' 're_sam_rate'...
    're_intaoxsp' 're_maxaoxsp' 're_tmaxaoxsp' 're_aotn' 're_aopinf' 're_aopn'...
    're_aofita' 're_aofitb' 're_aorsq' 're_prob' 're_version'...
    're_aitype' 're_hr' 're_sbp2' 're_intbapr' 're_maxbapr' 're_tmaxbapr'...
	're_intbaxsp' 're_maxbaxsp' 're_tmaxbap' 're_bafita' 're_bafitb'...
    're_barsq' 're_bapinf' 're_bapn' 're_ao_dpdt' 're_ba_dpdt'...
    're_pb_pf' 're_ri' 're_wf1i'  're_wf1t' 're_wf1a' 're_wbi' ...
    're_wbt' 're_wba' 're_wf2i'  're_wf2t' 're_wf2a'  're_wri' 're_rhoc' ...
    're_aosevr' 're_basevr' 're_quality'}; % header

% writetable
Results_table=cell2table(proc_var, 'VariableNames',header);
writetable(Results_table, xlsfile);

end
%% Tidy up
% rather than clearing the workspace for now I've left all non-redundant
% variables so that they can be used for debugging
clear p ptxt del dt i avlbeat b_qc calss_p xx yy zx sstxt newend extra ans a aa sstxt cuxwia;

% message at end
msg=string(no_of_files);
msg=strcat(msg, " files processed");
h1 = msgbox(msg, 'Done');
clear h1
