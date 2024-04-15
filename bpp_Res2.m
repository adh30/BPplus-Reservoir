%% batch analysis of BP+ data to perform pulse wave analysis, reservoir analysis and wave intensity analysis
%% Copyright 2019 Alun Hughes based on some original code by Kim Parker
% Also uses xml2struct.m by W. Falkena, ASTI, TUDelft, 21-08-2010 with additional
% modifications by A. Wanner, I Smirnov & Chao-Yuan Yeh and fill_between.m
% originally written by Ben Vincent, July 2014. Inspired by a function of
% the same name available in the Matplotlib Python library.

%%
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
% bpp_res2 (beta5) - in progress
% bpp_res2 (beta4) - modifications to ai_v1 based on suggestions from Richard Scott
% now renamed ai_v2.
% bpp_res2 (beta3) - beta version adapted to read old and new BPplus. Also
% incorporates additional information and suggestions from Richard Scott.
% Minor bug fixes to alternative folder setting, propagating default
% folder throughout the program and acknowledging use of 'fill_between.m'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Things to do
% - calculate other sphygmocor parameters (todo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% m files required to be in directory (in addition to bpp_Res2.m)
% ai_v2.m
% xml2struct.m
% BPfitres_v1.m
% kreservoir_v15.m
% fill_between.m
%% Constants
bRes_version='beta5';   % Version of bRes_bpp
kres_v='v15';           % Version tracking for reservoir fitting
headernumber=53;        % Headers for columns of results (see end)
mmHgPa = 133;           % P conversion for WIA
uconst=1;               % Empirical constant to convert normalized velocity to m/s
Npoly=3;                % Order of polynomial fit for sgolay
unprocessed_no = 0;     % Number of unprocessed files
Frame=9;                % Window length for sgolay based on (Rivolo et al.
% IEEE Engineering in Medicine and Biology Society
% Annual Conference 2014; 2014: 5056-9.

%% Select files
folder_name ='D:\BPPdata\'; % standard directory changed to reflect new xml files
% check that folder name exists and if not allows new folder to be chosen
if ~exist(folder_name, 'dir')
    answer = questdlg(folder_name + 'doesnt exist. Would you like to choose another folder?', ...
        'BPplus Data Folder','Yes', 'No [end]','Yes');
    % Handle response
    switch answer
        case 'Yes'
            folder_name = uigetdir;
            folder_name = strcat(folder_name,'\');  % adds \ to end of folder name for windows
        case 'No [end]'                             % end if no folder identified
            return
    end
end
file_lists=dir(fullfile(folder_name, '*.xml'));
no_of_files=length(file_lists);
% add an error trap here if no files in folder
if no_of_files==0
    f = errordlg('No data files to analyse in folder','File error');
    return
end
% set record number to 1 and extract filename
record_no=1;
% preallocate cell array
proc_var=cell(no_of_files,headernumber);

%% Start loop through all files
for file_number=1:no_of_files
    % refresh filename
    filename=file_lists(record_no).name;
    [data] = xml2struct([folder_name filename]);

    %file ID (i.e. filname without the .xml
    fileID = extractBefore( filename , '.xml');

    % progress bar
    % open waitbar
    msg = strcat('Processing file: ',filename);
    if record_no==1
        h = waitbar(record_no/no_of_files, msg);
    else
        waitbar(record_no/no_of_files, h, msg);
    end

    %% Read Data from CardioScope legacy xml
    if sum(strcmp(fieldnames(data), 'CardioScope')) == 1
        % extract CardioScope values
        ba_sbp=str2double(data.CardioScope.MeasDataLogger.Sys.Text);                % brachial systolic BP, mmHg
        dbp=str2double(data.CardioScope.MeasDataLogger.Dia.Text);                   % diastolic BP, mmHg
        ba_pp=ba_sbp-dbp;                                                           % brachial pulse pressure, mmHg
        map=str2double(data.CardioScope.MeasDataLogger.Mean.Text);                  % mean arterial pressure, mmHg
        hr=str2double(data.CardioScope.MeasDataLogger.Hr.Text);                     % heart rate, bpm
        samplerate=str2double(data.CardioScope.MeasDataLogger.SampleRate.Text);     % sample rate, Hz
        aosbp=str2double(data.CardioScope.Results.Result.aoSys.Text);               % cSBP calculated by BP+, mmHg
        aodbp=str2double(data.CardioScope.Results.Result.aoDia.Text);                % cDBP calculated by BP+, mmHg
        aopp=aosbp-aodbp;                                                           % cPP calculated by BP+, mmHg
        snr=str2double(data.CardioScope.Results.Result.SNR.Text);                   % Signal to noise ratio, dB
        ss_rmssd=str2double(data.CardioScope.Results.Result.RMSSD.Text);            % RMSSD from suprasystolic signal
        ss_ai=str2double(data.CardioScope.Results.Result.ssAI.Text);                % AI from suprasystolic signal
        ssdpdt=str2double(data.CardioScope.Results.Result.ssDpDtMax.Text);          % dp/dt from suprasystolic signal in uncorrected units
        ssHARM=str2double(data.CardioScope.Results.Result.ssHARM.Text);             % Normalized sAI
        ssPP=str2double(data.CardioScope.Results.Result.ssPP.Text);                 % Suprasystolic Pulse Pressure in the cuff, mmHg (PP in cuff is small)
        ssPPV=str2double(data.CardioScope.Results.Result.ssPPV.Text);               % Suprasystolic Pulse Pressure Variation, %
        ssRWTTFoot=str2double(data.CardioScope.Results.Result.ssRWTTFoot.Text);     % reflected wave transit time from foot of suprasystolic signal
        ssRWTTPeak=str2double(data.CardioScope.Results.Result.ssRWTTPeak.Text);     % reflected wave transit time from peak of suprasystolic signal
        ssSEP=str2double(data.CardioScope.Results.Result.ssSEP.Text);               % Systolic ejection period, s
        bppalgo=data.CardioScope.Results.Result.Attributes.algorithm_revision;      % Software algorithm
        Tx=split(data.CardioScope.Results.Result.ssAverageBeatPointsIdxs.Text,','); % Times of characteristic points?
        % Timings of peaks dont match waveforms - seems to be a bug or possibly
        % timings are from foot rather than beginning of waveform?
        % T1=str2double(Tn(2));                                                     % Time of 1st peak in samples
        % T2=str2double(Tn(3));                                                     % Time of inflection in samples
        % T3=str2double(Tn(4));                                                     % Time of 2nd peak in samples
        % T4=str2double(Tn(5));                                                     % Time of nadir of dichrotic notch
        datestring=data.CardioScope.MeasDataLogger.Attributes.datetime;             % Date as text string
        % brachial pulses
        ba_p_all=str2double(split(data.CardioScope.Results.Result.baEstimate.Text,','));
        % brachial average beat
        sstxt = split(data.CardioScope.Results.Result.ssAverageBeat.Text,',');
        % aortic pulses
        ao_p_all=str2double(split(data.CardioScope.Results.Result.aoEstimate.Text,','));
        % aortic average beat, ao_p_av
        a0=str2double(split(data.CardioScope.Results.Result.aoAverageBeat.Text,','));
        % create a double beat' to deal with the errors in definition of dbp
        a=[a0; a0];
        plot(a);
        % identify start and end of beat as minima
        [~, locmin1]=min(a);
        [~, locmin2]=min(a(locmin1+50:end));
        locmin2=locmin2+locmin1+50;
        plot(a(locmin1:locmin2));

        % filter derivative of new beat with SG to get rid of kinks at or around join
        aa = sgolayfilt(diff(a),Npoly,Frame); % filter the derivative - order and framelen defined above
        a1=cumsum(aa)-min(cumsum(aa))+min(a);    % reconstruct
        ao_p_av =a1(locmin1:locmin2-1)';         % crop to cycle
        clear a0 a aa a1;

        % start of pulses.
        baTransitTime = 0.18;
        if (contains(data.CardioScope.MeasDataLogger.Attributes.software_version,"038"))
            baTransitTime = 0.06;
        end
        aoOffset = round(baTransitTime * samplerate);
        ssBeatStartIdxs=str2double(split(data.CardioScope.Results.Result.ssBeatStartIdxs.Text,','));
        cPulseStartIndexes = ssBeatStartIdxs - aoOffset;
        for i = 1:length(cPulseStartIndexes)
            if cPulseStartIndexes(i)<0
                cPulseStartIndexes(i)=0;
            end
        end
    else
        %% Read Data from BPplus xml
        % extract BPplus values
        bppvers=data.BPplus.MeasDataLogger.Attributes.version;                      % Software version
        bppalgo=data.BPplus.Results.Result.Attributes.algorithm_revision;           % Software algorithm
        % Mode is only available in later versions of xml
        if isfield(data,'NibpMode')
            mode = string(data.BPplus.MeasDataLogger.NibpModeUsed.Text);            % Measurement mode
        end
        datestring=data.BPplus.MeasDataLogger.Attributes.datetime;                  % Date as text string
        ba_sbp=str2double(data.BPplus.MeasDataLogger.Sys.Text);                     % brachial systolic BP, mmHg
        dbp=str2double(data.BPplus.MeasDataLogger.Dia.Text);                        % diastolic BP, mmHg
        ba_pp=ba_sbp-dbp;                                                           % brachial pulse pressure, mmHg
        map=str2double(data.BPplus.MeasDataLogger.Map.Text);                        % mean arterial pressure, mmHg
        hr=str2double(data.BPplus.MeasDataLogger.Pr.Text);                          % heart rate, bpm
        samplerate=str2double(data.BPplus.MeasDataLogger.SampleRate.Text);          % sample rate, Hz
        aosbp=str2double(data.BPplus.Results.Result.cSys.Text);                     % cSBP calculated by BP+, mmHg
        aodbp=str2double(data.BPplus.Results.Result.cDia.Text);                     % cDBP calculated by BP+, mmHg
        aopp=aosbp-aodbp;                                                           % cPP calculated by BP+, mmHg
        snr=str2double(data.BPplus.Results.Result.SNR.Text);                        % Signal to noise ratio, dB
        ss_ai=str2double(data.BPplus.Results.Result.sAI.Text);                      % AI from suprasystolic signal
        ssdpdt=str2double(data.BPplus.Results.Result.sDpDtMax.Text);                % dp/dt from suprasystolic signal in uncorrected units
        ssPP=str2double(data.BPplus.Results.Result.sPP.Text);                       % Suprasystolic Pulse Pressure in the Cuff, mmHg (PP in cuff is small)
        ssPPV=str2double(data.BPplus.Results.Result.sPPV.Text);                     % Pulse pressure variation, %
        ssPRV=str2double(data.BPplus.Results.Result.sPRV.Text);                     % Pulse rate variation for selected pulses, ms (Root mean square of successive RR interval differences, RMSSD)
        ssRWTTFoot=str2double(data.BPplus.Results.Result.sRWTTFoot.Text);           % Reflected wave transit time from foot of suprasystolic signal
        ssRWTTPeak=str2double(data.BPplus.Results.Result.sRWTTPeak.Text);           % Reflected wave transit time from peak of suprasystolic signal
        ssSEP=str2double(data.BPplus.Results.Result.sSEP.Text)/1000;                % Systolic ejection period, s
        baTx=str2double(split(data.BPplus.Results.Result.sAveragePulsePointsIndexes.Text,',')); % Times of pulse points (1-6) see diagram
        ba_p_all=str2double(split(data.BPplus.Results.Result.baEstimate.Text,',')); % brachial pulses
        sstxt = split(data.BPplus.Results.Result.sAveragePulse.Text,',');           % brachial average beat ** not scaled ** [scaled later]
        ao_p_all=str2double(split(data.BPplus.Results.Result.cEstimate.Text,','));  % aortic pulses
        Selectedpulses = str2double(split(data.BPplus.Results.Result.sSelectedPulseIndexes.Text,','))+1; % selected pulses
        sPIndex = str2double(split(data.BPplus.Results.Result.sPulseStartIndexes.Text,','));

        %pPX brachial pulsatility index (PP/MAP)
        %cPX aortic pulsatility index (aoPP/MAP)

        % extract aortic data
        ao_p_av=str2double(split(data.BPplus.Results.Result.cAveragePulse.Text,','))';
        % start of pulses
        cPulseStartIndexes=str2double(split(data.BPplus.Results.Result.cPulseStartIndexes.Text,','));
    end
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
        quality='Unacceptable';
    end

    %% don't process poor or unacceptable files.
    if snr >=6 
        %% brachial pulses
        % replace values at start and end <DBP and >SBP with NaN
        % deal with low early values
        for i = 1:200
            if ba_p_all(i)<dbp
                ba_p_all(i) = NaN;
            end
        end
        % deal with high end values
        for i = length(ba_p_all)-200:length(ba_p_all)
            if ba_p_all(i)>ba_sbp
                ba_p_all(i) = NaN;
            end
        end
        % remove NaNs
        ba_p_all = ba_p_all(~isnan(ba_p_all));

        % calculate individual selected pulses for later display
        numgoodpulses = length(Selectedpulses);
        pulseindex = sPIndex-(sPIndex(1)-1);
        pulselengths = diff(pulseindex);
        pulsewaveforms=NaN(max(pulselengths),numgoodpulses);
        clear sPIndex

        for i = 1:numgoodpulses
            crop=ba_p_all(pulseindex(i):pulseindex(i+1));
            pulsewaveforms(1:length(crop),i)=crop;
        end

        % brachial average beat
        %sstxt = split(data.CardioScope.Results.Result.ssAverageBeat.Text,',');
        %b_avp_av=zeros(size(sstxt,1),size(sstxt,2));
        ba_p_av=str2double(sstxt);
        calss_p=ba_pp/(max(ba_p_av)-min(ba_p_av));
        ba_p_av=dbp+(ba_p_av*calss_p);
        ssdpdt=ssdpdt*calss_p;                                                      % correcting to mmHg/s

        % aortic pulses
        ao_p_all=ao_p_all(cPulseStartIndexes(1):cPulseStartIndexes(length(cPulseStartIndexes)));
        %%ao_p_all=str2double(split(data.CardioScope.Results.Result.aoEstimate.Text,','));
        % replace values at start and end <DBP and >SBP with NaN
        % deal with low early values
        for i = 1:200
            if ao_p_all(i)<dbp
                ao_p_all(i) = NaN;
            end
        end
        % deal with high end values
        for i = length(ao_p_all)-200:length(ao_p_all)
            if ao_p_all(i)>aosbp
                ao_p_all(i) = NaN;
            end
        end
        % remove NaNs
        ao_p_all = ao_p_all(~isnan(ao_p_all));

        % calculate aoAIx, Pi, Tfoot, Ti (T1), Tpeak(T2)
        [ao_ai, aoPi, aoTfoot, aoTi, aoTmax, aoTypetxt] = ai_v2(ao_p_av, samplerate);
        ao_Tr = aoTi-aoTfoot;
        aoPidbp = aoPi-dbp;

        %% calculate aortic dp/dt
        ao_dpdt=max(diff(ao_p_av))*samplerate;

        %% do reservoir calculations for aortic pressure
        [aoTn_av, aoPinf_av, aoP_av,aoPr_av,aoPn_av, aofita_av,...
            aofitb_av, aorsq_av]=BPfitres_v1(ao_p_av,samplerate);
        aoPxs=aoP_av-aoPr_av;

        %% do reservoir calculations for brachial pressure
        % assign p to pass to function
        [baTn_av, baPinf_av, baP_av,baPr_av,baPn_av, bafita_av,...
            bafitb_av, barsq_av]=BPfitres_v1(ba_p_av',samplerate);

        % calculate Pxs = P - Pres
        baPxs=baP_av-baPr_av;

        % pAI
        [~, Tdiffpeaks]=findpeaks(diff(ba_p_av), 'Npeaks', 3);    % peaks of dP
        ba_p2=ba_p_av(Tdiffpeaks(2));
        ba_ai=100*(ba_p2-dbp)/(ba_sbp-dbp); % peripheral augmentation index (%) uses definition in Munir S et al. Hypertension 2008; 51(1): 112-8.
        ba_dpdt=max(diff(ba_p_av))*samplerate;

        %% calculate SEVR (aortic)
        % aoTn and baTn are not necessarily the same - which is better is uncertain
        % but on the basis that it's aortic not brachial which is the cardiac load
        % we'll use that.
        % sub endocardial viability ratio (SEVR)
        % SEVR = diastolic pressure-time integral(DPTI)/systolic pressure-time integral(SPTI)
        lsys=round(aoTn_av*samplerate);
        ao_spti=trapz(ao_p_av(1:lsys));
        ao_dpti=trapz(ao_p_av(lsys+1:end));
        ao_sevr=ao_dpti/ao_spti;                    % SEVR (aortic)
        kliuinverse=ao_dpti/(ao_dpti+ao_spti);      % inverse of Liu constant, k
        pmean_sys=ao_spti/lsys;                     % mean pressure of systole
        pmean_dia=ao_dpti/(length(ao_p_av)-lsys);   % mean pressure of diastole

        %% Create estimates of Pf and Pb using the assumption that Pb = cPr/2
        aoPb_av=(aoPr_av-min(aoP_av))/2;
        aoPf_av=aoP_av-aoPb_av-min(aoP_av);
        Pb_Pf=max(aoPb_av)/max(aoPf_av);
        RI=max(aoPb_av)/(max(aoPf_av)+ max(aoPb_av));

        %% wave intensity analysis (only done on central P)
        % convert Pxs to flow velocity and assume peak velocity, U = 1m/s.
        % Estimate of U based on Lindroos, M., et al. J Am Coll Cardiol
        % 1993; 21: 1220-1225.
        cu=uconst*(aoPxs/max(aoPxs));
        u=cu;
        u(9:end)=cu(1:end-8);
        u(1:3)=0;
        u(4:9)= nan;
        u=fillmissing(u,'spline')';
        u=u';   %
        % calculate derivatives
        [b,g]=sgolay(Npoly,Frame);   % Calculate S-G coefficients
        HalfWin=((Frame+1)/2) -1;
        N=length(aoP_av);
        dp=zeros(N,1); duxs=dp;
        for n=(Frame+1)/2:N-(Frame+1)/2
            % 1st differential
            dp(n)=dot(g(:,2),aoP_av(n-HalfWin:n+HalfWin));  % pressure difference
            duxs(n)=dot(g(:,2),u(n-HalfWin:n+HalfWin));     % velocity difference
        end

        di=dp.*duxs;
        di=di*mmHgPa*length(dp)^2;      % units fixed - now in W/m2 per cycle^2

        % new peak detection algorithm (version beta 5)
        minpeak=max(di)/20;             % changed to 20 - entirely arbitrary but represents 5%.
        [~,lsys]=min(dp);               % restrict analysis to systole
        %lsys=round(lsys)+5;            % round and add 5 samples to give a margin for error for duration of systole
        lsys=ceil(lsys);                % roundup duration of systole

        % Wf1 defined as the first (substantial) positive wave intensity peak
        % to avoid occasional odd double peaks for Wf1 use maximum peak to guide
        % search for first peak
        [dippks(1),diplocs(1), dipw(1)]=findpeaks(di(1:lsys+10), 'NPeaks',1,'MinPeakHeight',max(di)*.9);

        % Wb defined as the largest negative peak that follows Wf1
        [dimpks,dimlocs,dimw]=findpeaks(-di(1:lsys), 'NPeaks',1,'MinPeakHeight',0.7*max(-di)); % find one dI- peaks (Wb)

        % Wf2 defined as largest positive peak when dp is negative
        [dippks(2),diplocs(2), dipw(2)]=findpeaks(di(lsys-20:lsys+20), 'NPeaks',1, 'SortStr','descend'); % find 2nd dI+ peaks (Wf2) on the asssumption it follows Wf2 and allowing 5 samples beyond length of systole.
        diplocs(2)=diplocs(2)+(lsys-20); % add lsys-20 since we are only searching from lsys-20

        % calculate peak time (s)
        dipt=diplocs/samplerate;
        dimt=dimlocs/samplerate;

        % calculate areas
        % For a Gaussian curve (assumed) the area is 1.06447*height*width
        diparea=1.06447*dippks.*dipw;
        dimarea=1.06447*dimpks.*dimw;
        wri=dimarea/diparea(1);

        % % error trap when W2 is unmeasureable
        %  if length(dippks)==1
        %     dippks(2)=0;
        %     dipt(2)=0;
        %     diparea(2)=0;
        %  end

        % Estimate c (wavespeed) as k*dP/du where k is empirical constant
        rhoc=max(aoPxs)*mmHgPa/1000; % units (m/s)

        %% possible problems with fits
        % ****** MORE TO DO HERE *********************
        if aoPinf_av>=dbp || aoPinf_av<-12 || aofita_av<=0 || aofitb_av<=0 || aorsq_av<0.9
            prob=1;
        else
            prob =0;
        end
        % *************should write error salvage sometime!

        %% make figures and data subfolders
        figfolder=strcat(folder_name, 'figures\');
        datafolder=strcat(folder_name,'results\');
        if ~exist(figfolder, 'dir')
            mkdir(figfolder);
        end
        if ~exist(datafolder, 'dir')
            mkdir(datafolder);
        end
        % save separate jpg (for viewing)
        % wmf deprecated due to security issues, svg not functional
        %expression = 'xml';
        %replace2 ='jpg';
        %jpgfile = regexprep(filename,expression,replace2);

        %% Print figures
        % make a time variable for printing
        Time_cycles=(1:length(pulsewaveforms))/samplerate;
        % Aortic pressure and reservoir
        figure('Name',fileID,'visible','on');
        subplot(2,2,1);
        plot(Time_cycles,pulsewaveforms);
        xlabel('Time (s)')
        ylabel('BP (mmHg)')
        title('Pulse traces')
        box off;
        % Reservoir
        subplot(2,2,2);
        % Time_av=(1:length(aoPf_av))/samplerate;
        % plot(Time_av, aoPf_av,'b', Time_av, aoPb_av,'r', Time_av, aoP_av-min(aoP_av),'k');
        % xlabel('Time (s)')
        % ylabel('BP (mmHg)')
        % title('aortic P, Pf and Pb')
        Time_av=(1:length(aoPf_av))/samplerate;
        plot(Time_av, aoPr_av-min(aoP_av),'b', Time_av, aoPxs,'r', Time_av, aoP_av-min(aoP_av),'k');
        xlabel('Time (s)')
        ylabel('BP (mmHg)')
        title('aortic P, Pres and Pxs')
        % SEVR
        subplot(2,2,3)
        x = 1:length(ao_p_av);
        % xScale=5;
        plot(x,ao_p_av); hold on;
        xlabel('Samples')
        ylabel('BP (mmHg)')
        base=min(ao_p_av);
        dicroticNotchIdx=fix(ssSEP*samplerate);         % This uses the BP+ determined time of dicrotic notch.
        if (lsys <= dicroticNotchIdx)
            fill_between(1:length(ao_p_av),ao_p_av, base, x <= lsys, 'FaceColor', [0.8 0.8 1]);
            fill_between(1:length(ao_p_av),ao_p_av, base, x >= lsys & x <= dicroticNotchIdx , 'FaceColor', [1 0.8 1]);
            fill_between(1:length(ao_p_av),ao_p_av, base, x >= dicroticNotchIdx, 'FaceColor', [1 0.8 0.8]);
        else
            fill_between(1:length(ao_p_av),ao_p_av, base, x <= dicroticNotchIdx, 'FaceColor', [0.8 0.8 1]);
            fill_between(1:length(ao_p_av),ao_p_av, base, x >= dicroticNotchIdx & x <= lsys , 'FaceColor', [1 1 0.8]);
            fill_between(1:length(ao_p_av),ao_p_av, base, x >= lsys, 'FaceColor', [1 0.8 0.8]);
        end
        x = [lsys,lsys];
        y = [base,ao_p_av(lsys)];
        plot(x,y);
        x = [dicroticNotchIdx,dicroticNotchIdx];
        y = [base,ao_p_av(dicroticNotchIdx)];
        plot(x,y);
        %sevr_title = sprintf('%s SEVR=%3.0f %% cAIx=%3.0f %% cAP=%3.0f %% cPP=%3.0f %%',filename,ao_sevr*100,ao_ai,aosbp-aoPi,aopp);
        title('SEVR')
        %cAIx_txt = sprintf('\\leftarrow T1',aoTmax);
        cAIx_txt = '\leftarrow T1';
        text((aoTi*samplerate)+1,ao_p_av(aoTi*samplerate),cAIx_txt,'FontSize', 6);
        % if to catch when lsys==dicroticNotchIdx so they dont overwrite.
        if lsys~=dicroticNotchIdx
            maxnegdPdt_txt = '\leftarrow max(-dP/dt)';
            text((lsys+1),ao_p_av(lsys),maxnegdPdt_txt,'FontSize', 6);
        else
            maxnegdPdt_txt = 'max(-dP/dt) =';
            text((lsys+1),ao_p_av(lsys)+.075*(max(ao_p_av)-base),maxnegdPdt_txt,'FontSize', 6);
        end
        dicroticNotchIdx_txt = '\leftarrow notch';
        text(dicroticNotchIdx+1,ao_p_av(dicroticNotchIdx),dicroticNotchIdx_txt,'FontSize', 6);
        text((lsys/2),base+10,'SPTI','FontSize', 6, HorizontalAlignment='center');
        text((lsys+lsys/2),base+10,'DPTI','FontSize', 6);
        % sevr_txt = sprintf('SEVR=%3.0f %%',ao_sevr*100);
        % text((lsys+lsys),ao_p_av(lsys)-10,sevr_txt,'FontSize', 6);

        % WI
        subplot(2,2,4)
        TimeDI=(1:length(di))/samplerate;
        plot (TimeDI, di); hold on;                         % to allow new length
        plot(dipt(1),dippks(1),'ko');
        plot(dimt,-dimpks,'ro');
        plot(dipt(2),dippks(2),'ks');
        xlabel('Time (s)')
        ylabel('dI (W/m^2/cycle^2)')
        title('Wave intensity, dI')

        set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',8);
        % Save
        jpgfile = regexprep(filename,'.xml','figs.jpg');
        print ('-djpeg', '-r300' , [figfolder jpgfile]);

        % clear and close figures
        clear f1 f2
        figs =  findobj('type','figure');
        close(figs);
        clear figs;

        %% Calculate and save additional results
        [~,max_tpb]= max(ba_p_av);          % maximum forward pressure
        [~,max_tpa]= max(aoP_av);           % maximum backward pressure
        [maxPra,max_tra]= max(aoPr_av);     % maximum aortic reservoir pressure and time
        [maxPxsa,max_txsa]= max(aoPxs);     % maximum aortic excess pressure and time
        [maxPrb,max_trb]= max(baPr_av);     % maximum brachial reservoir pressure and time
        [maxPxsb,max_txsb]= max(baPxs);     % maximum brachial excess pressure and time

        % wasted LV pressure energy (Ew) = 2.09Δtr (Ps − Pi). Nichols WW.
        % Clinical measurement of arterial stiffness obtained from
        % noninvasive pressure waveforms. Am J Hypertens 2005; 18(1 Pt 2): 3S-10S.



        %% write variables
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
        proc_var{record_no,21}=kres_v;                                          % version kreservoir
        proc_var{record_no,22}=aoTypetxt;                                       % AI type
        proc_var{record_no,23}=hr;                                              % Heart rate, bpm
        proc_var{record_no,24}=ba_p2;                                         % SBP2, mmHg
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
        proc_var{record_no,52}='bpp_Res2 (beta2)';                              % version of bpp_Res
        proc_var{record_no,53}=quality;                                         % quality index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % Error trap for low SNR and message
        unprocessed_no = unprocessed_no+1;
        disp(strcat(filename, ' not processed due to poor quality'));
    end

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
    % 52 'empty'
    % 53 'quality'
    % % % %
    xlsfile=strcat(folder_name, 'results\resdata.xls');
    header = {'re_file' 're_basbp' 're_tbasbp' 're_minp' 're_intaopr' 're_maxaopr'...
        're_tmaxaopr' 're_intaoprlessdbp' 'date' 're_sam_rate'...
        're_intaoxsp' 're_maxaoxsp' 're_tmaxaoxsp' 're_aotn' 're_aopinf' 're_aopn'...
        're_aofita' 're_aofitb' 're_aorsq' 're_prob' 're_kres'...
        're_aitype' 're_hr' 're_sbp2' 're_intbapr' 're_maxbapr' 're_tmaxbapr'...
        're_intbaxsp' 're_maxbaxsp' 're_tmaxbap' 're_bafita' 're_bafitb'...
        're_barsq' 're_bapinf' 're_bapn' 're_ao_dpdt' 're_ba_dpdt'...
        're_pb_pf' 're_ri' 're_wf1i'  're_wf1t' 're_wf1a' 're_wbi' ...
        're_wbt' 're_wba' 're_wf2i'  're_wf2t' 're_wf2a'  're_wri' 're_rhoc' ...
        're_aosevr' 're_version' 're_quality'}; % header

    % % writetable
    Results_table=cell2table(proc_var, 'VariableNames',header);
    writetable(Results_table, xlsfile);
end
%% Tidy up
% rather than clearing the workspace for now I've left all non-redundant
% variables so that they can be used for debugging
% clear cuxwia;

% message at end
msg=string(no_of_files);
msg=strcat(msg, " files processed");
msgbox(msg, 'Done');


