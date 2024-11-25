%% batch analysis of BP+ data to perform pulse wave analysis, reservoir analysis and wave intensity analysis
%% Copyright 2019 Alun Hughes based on some original code by Kim Parker
% Also uses xml2struct.m by W. Falkena, ASTI, TUDelft, 21-08-2010 with additional
% modifications by A. Wanner, I Smirnov & Chao-Yuan Yeh
% and
% fill_between.m
% originally written by Ben Vincent, July 2014. Inspired by a function of
% the same name available in the Matplotlib Python library.
% and
% MMINTERP by D.C. Hanselman, Jan 1994

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
% bpp_res2 (beta5) - Expand calculated variables so they match key Sphygmocor variables; 
% - Use separate functions to read Cardioscope and BP+ files;
% - Improve detection of end of systole; 
% - Fix bug in Wf1 and Wf2 peak identification 
% - Better identification of problems with fits 
% - Some error salvage
% - Fix bugs in figures
% - Cut number of figures to one per BP+ recording

% bpp_res2 (beta4) - modifications to ai_v1 based on suggestions from
% Richard Scott now renamed ai_v2.
% bpp_res2 (beta3) - beta version adapted to read old and new BPplus. Also
% incorporates additional information and suggestions from Richard Scott.
% Minor bug fixes to alternative folder setting, propagating default
% folder throughout the program and acknowledging use of 'fill_between.m'
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
headernumber=54;        % Headers for columns of results (see end)
mmHgPa = 133.322;           % P conversion for WIA
uconst=1;               % Empirical constant to convert normalized velocity to m/s
unprocessed_no = 0;     % Number of unprocessed files
Npoly=3;                % Order of polynomial fit for sgolay
Frame=9;                % Window length for sgolay based on (Rivolo et al.
% IEEE Engineering in Medicine and Biology Society
% Annual Conference 2014; 2014: 5056-9).
%% Select files and folder
% folder from json
jtext = fileread('bppconfig.json');
jdata = jsondecode(jtext);
%file_lists=dir(fullfile(string(data.folder_name), '*.*'));
folder_name = jdata.folder_name;
disp(['Source folder = ' folder_name]);

%folder_name ='D:\BPPdata\'; % standard directory changed to reflect new xml files
% check that folder name exists and if not allows new folder to be chosen
if ~exist(folder_name, 'dir')
    answer = questdlg([folder_name ' doesnt exist. Would you like to choose another folder?'], ...
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
% ensure subdirectories exist to be written to
% TBD

% identify xml files
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

    % progress bar
    msg = strcat(['Processing file ',num2str(record_no), ' of ', num2str(no_of_files)]);
    if record_no==1
        h = waitbar(record_no/no_of_files, msg);
    else
        waitbar(record_no/no_of_files, h, msg);
    end

    %% Read BP+ xml file (Both legacy 'CardioScope' and current 'BPplus' xml are supported.)
    [data, metadata, ss, ba, ao]=read_BPplus(folder_name, filename, Npoly, Frame);
    samplerate=metadata.samplerate;                                 % assign to local variable to simplify use of this common value

    %% don't process poor or unacceptable files.
    if metadata.snr >=6
        %% process acceptable files
        % quality of BP tracing
        if metadata.snr < 9
                quality = 'acceptable';
        elseif metadata.snr >=9 && metadata.snr<12
                quality = 'good';
        else
                quality = 'excellent';
        end

        % fix erroneous zeros in pulsewaveforms
        ss.pulsewaveforms(ss.pulsewaveforms == 0) = NaN;
        % calculate aoAIx, Pi, Tfoot, Ti, T1, Tmax, T2
        % T1 and T2 calculated according to Qasem & Avolio 2008 which may
        % differ from estimate of Ti by Kelly

        % Aortic waveform features including AIx
        [ao_ai, ao_Pi, ao_Tfoot, ao_Ti, ao_Tmax, ao_Typetxt] = ai_v2(ao.p_av, samplerate);
        ao_Tr = ao_Ti-ao_Tfoot;
        ao_Piba_dbp = ao_Pi-ba.dbp;
        if ao_Typetxt=="Type A"
            ao_t1= ao_Ti;
            ao_t2 = ao_Tmax;
        else
            ao_t1=ao_Tmax;
            ao_t2=ao_Ti;
        end
        ao_p1=ao.p_av(round(ao_t1*samplerate));
        ao_p2=ao.p_av(round(ao_t2*samplerate));
        ao_ap=ao_p2-ao_p1;

        %% calculate aortic dp/dt [mmHg/s]
        ao_dpdt=max(diff(ao.p_av))*samplerate;

        %% do reservoir calculations for aortic pressure
        [aoTn_av, aoPinf_av, aoP_av, aoPr_av, aoPn_av, aofita_av, aofitb_av, aorsq_av]=BPfitres_v1(ao.p_av, samplerate);
        aoPxs=aoP_av-aoPr_av;

        %% do reservoir calculations for brachial pressure
        [baTn_av, baPinf_av,baP_av, baPr_av,baPn_av, bafita_av, bafitb_av, barsq_av]=BPfitres_v1(ba.p_av',metadata.samplerate);
        % calculate Pxs = P - Pres where P is the processed brachial
        % pressure which has the same length as the reservoir pressure
        baPxs=baP_av-baPr_av;

        % pAI and fiducial points on brachial waveform
        [~, Tdiffpeaks]=findpeaks(diff(baP_av), 'Npeaks', 3);               % peaks of dP
        ba_p1=baP_av(Tdiffpeaks(1));
        ba_t1=Tdiffpeaks(1)/samplerate;
        ba_p2=baP_av(Tdiffpeaks(2));
        ba_t2=Tdiffpeaks(2)/samplerate;
        ba_ai=100*(ba_p2-ba.dbp)/(ba.sbp-ba.dbp);                           % peripheral augmentation index (%) uses definition in Munir S et al. Hypertension 2008; 51(1): 112-8.
        ba_dpdt=max(diff(baP_av))*samplerate;                               % mmHg/s
        [~,ba_tmax]=max(baP_av);
        ba_tmax=ba_tmax/samplerate;                                         %s

        %% calculate SEVR (aortic)
        % aoTn and baTn are not necessarily the same - which is better is uncertain
        % but on the basis that it's aortic not brachial which is the cardiac load
        % we'll use that.
        % sub endocardial viability ratio (SEVR)
        % SEVR = diastolic pressure-time integral(DPTI)/systolic pressure-time integral(SPTI)
        lsys=round(aoTn_av*samplerate);
        ao_spti=trapz(ao.p_av(1:lsys));
        ao_dpti=trapz(ao.p_av(lsys+1:end));
        ao_sevr=ao_dpti/ao_spti;                    % SEVR (aortic)
        kliuinverse=ao_dpti/(ao_dpti+ao_spti);      % inverse of Liu constant, k
        ao_pmean_sys=ao_spti/lsys;                     % mean pressure of systole
        ao_pmean_dia=ao_dpti/(length(ao.p_av)-lsys);   % mean pressure of diastole

        %% Create estimates of Pf and Pb using the assumption that Pb = cPr/2
        aoPb_av=(aoPr_av-min(aoP_av))/2;
        aoPf_av=aoP_av-aoPb_av-min(aoP_av);
        [aoPbmax,aot_Pbmax] = max(aoPb_av);         % maximum backward pressure, mmHg & samples
        aot_Pbmax=aot_Pbmax/samplerate;             % time maxumum backward pressure, s
        [aoPfmax, aot_Pfmax] = max(aoPf_av);        % maximum forward pressure, mmHg & samples
        aot_Pfmax=aot_Pfmax/samplerate;             % time maxumum forward pressure, s
        aoPb_Pf=aoPbmax/aoPfmax;                    % Pb/Pf
        ao_ri=aoPbmax/(aoPfmax + aoPbmax);          % Reflection index

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
        if aoPinf_av>=ba.dbp || aoPinf_av<-12 || aofita_av<=0 || aofitb_av<=0 || aorsq_av<0.9
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
        Time_cycles=(1:length(ss.pulsewaveforms))/samplerate;
        % Aortic pressure and reservoir
        %figure('Name',metadata.fileID,'visible','on');
        figure('Name',metadata.fileID,'visible','off');
        subplot(2,2,1);
        plot(Time_cycles,ss.pulsewaveforms);
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
        x = 1:length(ao.p_av);
        % xScale=5;
        plot(x,ao.p_av); hold on;
        xlabel('Samples')
        ylabel('BP (mmHg)')
        base=min(ao.p_av);
        dicroticNotchIdx=fix(ss.sep*samplerate);         % This uses the BP+ determined time of dicrotic notch (rounded down)
        if (lsys <= dicroticNotchIdx)
            fill_between(1:length(ao.p_av),ao.p_av, base, x <= lsys, 'FaceColor', [0.8 0.8 1]);
            fill_between(1:length(ao.p_av),ao.p_av, base, x >= lsys & x <= dicroticNotchIdx , 'FaceColor', [1 0.8 1]);
            fill_between(1:length(ao.p_av),ao.p_av, base, x >= dicroticNotchIdx, 'FaceColor', [1 0.8 0.8]);
        else
            fill_between(1:length(ao.p_av),ao.p_av, base, x <= dicroticNotchIdx, 'FaceColor', [0.8 0.8 1]);
            fill_between(1:length(ao.p_av),ao.p_av, base, x >= dicroticNotchIdx & x <= lsys , 'FaceColor', [1 1 0.8]);
            fill_between(1:length(ao.p_av),ao.p_av, base, x >= lsys, 'FaceColor', [1 0.8 0.8]);
        end
        x = [lsys,lsys];
        y = [base,ao.p_av(lsys)];
        plot(x,y);
        x = [dicroticNotchIdx,dicroticNotchIdx];
        y = [base,ao.p_av(dicroticNotchIdx)];
        plot(x,y);
        %sevr_title = sprintf('%s SEVR=%3.0f %% cAIx=%3.0f %% cAP=%3.0f %% cPP=%3.0f %%',filename,ao_sevr*100,ao_ai,ao.sbp-aoPi,ao.pp);
        title('SEVR')
        %cAIx_txt = sprintf('\\leftarrow T1',aoTmax);
        cAIx_txt = '\leftarrow T1';
        text((ao_Ti*samplerate)+1,ao.p_av(ao_Ti*samplerate),cAIx_txt,'FontSize', 6);
        % if to catch when lsys==dicroticNotchIdx so they dont overwrite.
        if lsys~=dicroticNotchIdx
            maxnegdPdt_txt = '\leftarrow max(-dP/dt)';
            text((lsys+1),ao.p_av(lsys),maxnegdPdt_txt,'FontSize', 6);
        else
            maxnegdPdt_txt = 'max(-dP/dt) =';
            text((lsys+1),ao.p_av(lsys)+.075*(max(ao.p_av)-base),maxnegdPdt_txt,'FontSize', 6);
        end
        dicroticNotchIdx_txt = '\leftarrow notch';
        text(dicroticNotchIdx+1,ao.p_av(dicroticNotchIdx),dicroticNotchIdx_txt,'FontSize', 6);
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
        pp_ar=ba.pp/ao.pp;                                                      % pulse pressure amplification ratio
        [maxPra,max_tra] = max(aoPr_av);                                        % maximum aortic reservoir pressure and time, mmHg & samples
        max_tra=max_tra/samplerate;                                             % s
        [maxPxsa,max_txsa] = max(aoPxs);                                        % maximum aortic excess pressure and time, mmHg & samples
        max_txsa=max_txsa/samplerate;                                           % s
        [maxPrb,max_trb] = max(baPr_av);                                        % maximum brachial reservoir pressure and time, mmHg & samples
        max_trb=max_trb/samplerate;                                             % s
        [maxPxsb,max_txsb] = max(baPxs);                                        % maximum brachial excess pressure and time, mmHg & s
        max_txsb=max_txsb/samplerate;                                           % s
        ao_ai75=ao_ai+.481*ba.hr-36.1;                                          % AIx heart rate adjusted to heart rate 75
        ao_dd=60/ba.hr-aoPn_av/1000;                                            % central diastolic duration, s
        e_w = pi/2*(aoTn_av-ao_Tr)*(ao.sbp-ao_Pi)*mmHgPa;                       % wasted LV pressure energy (Ew) = pi/2*Δtr (Ps − Pi)* conversion from mmHg to Pa. 
                                                                                % Nichols WW. Clinical measurement of arterial stiffness obtained 
                                                                                % from noninvasive pressure waveforms. Am J Hypertens 2005; 
                                                                                % 18(1 Pt 2): 3S-10S.

        %% write variables
        proc_var{record_no,1}=filename;                                         % filename
        proc_var{record_no,2}=metadata.datestring;                              % Date as text string
        proc_var{record_no,3}=metadata.bppvers;                                 % BP+ software version
        proc_var{record_no,4}=metadata.bppalgo;                                 % BP+ algorithm
        proc_var{record_no,5}=bRes_version;                                     % bRes version
        proc_var{record_no,6}=kres_v;                                           % kreservoir version
        proc_var{record_no,7}=samplerate;                                       % sampling rate, 1/sec
        proc_var{record_no,8}=ba.sbp;                                           % brachial SBP, mmHg
        proc_var{record_no,9}=ba_tmax;                                          % time max brachial P (baSBP), s
        proc_var{record_no,10}=ba.dbp;                                          % min P, ba.dbp, mmHg
        proc_var{record_no,11}=ba.hr;                                           % Heart rate, bpm
        proc_var{record_no,12}=ba.map;                                          % Mean arterial pressure, mmHg
        proc_var{record_no,13}=ba.pp;                                           % Brachial pulse pressure, mmHg
        proc_var{record_no,14}=ao.sbp;                                          % Aortic systolic pressure, mmHg
        proc_var{record_no,15}=ao.dbp;                                          % Aortic diastolic pressure, mmHg
        proc_var{record_no,16}=ao.pp;                                           % Aortic pulse pressure, mmHg
        proc_var{record_no,17}=metadata.snr;                                    % Signal to Noise Ratio
        proc_var{record_no,18}=ss.prv;                                          % Pulse rate variability from suprasystolic signal, ms
        proc_var{record_no,19}=ss.ai;                                           % Augmentation index from suprasystolic signal, ms
        proc_var{record_no,20}=ss.ppv;                                          % Pulse pressure variation from suprasystolic signal, mmHg
        proc_var{record_no,21}=ss.RWTTFoot;                                     % Reflected wave transit time from foot of suprasystolic signal
        proc_var{record_no,22}=ss.RWTTPeak;                                     % Reflected wave transit time from peak of suprasystolic signal
        proc_var{record_no,23}=ss.sep;                                          % Systolic ejection period from suprasystolic signal (notch)
        proc_var{record_no,24}=quality;                                         % quality of trace based on SNR
        proc_var{record_no,25}=ba_t1;                                           % brachial T1, s
        proc_var{record_no,26}=ba_p1;                                           % brachial P1, mmHg
        proc_var{record_no,27}=ba_t2;                                           % brachial T2, s
        proc_var{record_no,28}=ba_p2;                                           % brachial P2 (SBP2), mmHg
        proc_var{record_no,29}=ba_ai;                                           % Peripheral augmentation index (P2/P1)
        proc_var{record_no,30}=baPn_av;                                        % brachial end systolic pressure, mmHg
        proc_var{record_no,31}=ba_dpdt;                                         % brachial dp/dt, mmHg/s
        proc_var{record_no,32}=aoTn_av;                                         % aortic ejection duration, s
        proc_var{record_no,33}=aoPn_av;                                         % aortic end-systolic pressure, mmHg/s
        proc_var{record_no,34}=ao_p1;                                           % aortic pressure at T1, mmHg
        proc_var{record_no,35}=ao_p2;                                           % aortic pressure at T2, mmHg
        proc_var{record_no,36}=ao_Tr;                                           % aortic pulse reflection time, s
        proc_var{record_no,37}=ao_Typetxt;                                      % AI type
        proc_var{record_no,38}=pp_ar;                                           % Pulse pressure amplification ratio (brachial/central)
        proc_var{record_no,39}=ao_ap;                                           % aortic augmented pressure, mmHg
        proc_var{record_no,40}=ao_pmean_sys;                                    % aortic mean pressure of systole, mmHg
        proc_var{record_no,41}=ao_pmean_dia;                                    % aortic mean pressure of diastole, mmHg
        proc_var{record_no,42}=ao_spti;                                         % aortic tension time index (spti), mmHg.s
        proc_var{record_no,43}=ao_dpti;                                         % aortic dpti, mmHg.s
        proc_var{record_no,44}=ao_sevr;                                         % aortic Sub-Endocardial Viability Ratio (SEVR, Buckberg Index)
        proc_var{record_no,45}=ao_dd;                                           % aortic diastolic duration, s
        proc_var{record_no,46}=ao_ai;                                           % aortic augmentation index (Central Aug/PH), % 
        proc_var{record_no,47}=ao_ai75;                                         % aortic heart rate corrected augmentation index, %
        proc_var{record_no,48}=ao_Ti;                                           % aortic time of shoulder/inflection point, s
        proc_var{record_no,49}=e_w;                                             % aortic wasted work, Pa.s
        proc_var{record_no,50}=ao_dpdt;                                         % aortic dp/dt, mmHg/s
        proc_var{record_no,51}=aoPbmax;                                         % aortic maximum backwards pressure, mmHg
        proc_var{record_no,52}=aot_Pbmax;                                       % aortic time maximum backwards pressure, s
        proc_var{record_no,53}=aoPfmax;                                         % aortic maximum forwards pressure, mmHg
        proc_var{record_no,54}=aot_Pfmax;                                       % aortic time maximum forward pressure, s
        proc_var{record_no,55}=aoPb_Pf;                                         % aortic Pb/Pf (ratio of backward to forward pressure - also known as reflection magnitude. RM)
        proc_var{record_no,56}=ao_ri;                                           % aortic reflection index (Pb/(Pb+Pf))
        proc_var{record_no,57}=sum(aoPr_av)/samplerate;                         % integral aoPres, mmHg.s
        proc_var{record_no,58}=maxPra;                                          % max aortic Pres, mmHg
        proc_var{record_no,59}=max_tra/samplerate;                              % Time max aortic Pres, s                  ****************RENAME
        proc_var{record_no,60}=sum(aoPr_av-ba.dbp)/samplerate;                  % integral aortic Pres-diastolic, mmHg.s
        proc_var{record_no,61}=sum(aoPxs)/samplerate;                           % Integral excess pressure, mmHg.s
        proc_var{record_no,62}=maxPxsa;                                         % max excess P, mmHg
        proc_var{record_no,63}=max_txsa/samplerate;                             % time of max excess P, sec
        proc_var{record_no,64}=aoPinf_av;                                       % aortic P infinity, mmHg
        proc_var{record_no,65}=aofita_av;                                       % aortic rate constant A (ka), 1/sec
        proc_var{record_no,66}=aofitb_av;                                       % aortic rate constant B (kb), 1/sec
        proc_var{record_no,67}=aorsq_av;                                        % aortic R2 for diastolic fit
        proc_var{record_no,68}=prob;                                            % likely problem with fit
        proc_var{record_no,69}=sum(baPr_av)/samplerate;                         % integral baPres, mmHg.s
        proc_var{record_no,70}=maxPrb;                                          % max baPres, mmHg
        proc_var{record_no,71}=max_trb/samplerate;                              % Time max Pres, sec
        proc_var{record_no,72}=sum(baPxs)/samplerate;                           % Integral ba excess pressure, mmHg.s
        proc_var{record_no,73}=maxPxsb;                                         % max ba excess P, mmHg
        proc_var{record_no,74}=max_txsb/samplerate;                             % time of max ba excess P, sec
        proc_var{record_no,75}=bafita_av;                                       % brachial rate constant A (ka), 1/sec
        proc_var{record_no,76}=bafitb_av;                                       % brachial rate constant B (kb), 1/sec
        proc_var{record_no,77}=barsq_av;                                        % brachial R2 for diastolic fit
        proc_var{record_no,78}=baPinf_av;                                       % brachial Pinf, mmHg
        proc_var{record_no,79}=dippks(1);                                       % W1 intensity
        proc_var{record_no,80}=dipt(1);                                         % W1 time
        proc_var{record_no,81}=diparea(1);                                      % W1 area
        proc_var{record_no,82}=dimpks;                                          % W-1 intensity
        proc_var{record_no,83}=dimt;                                            % W-1 time
        proc_var{record_no,84}=dimarea;                                         % W-1 area
        proc_var{record_no,85}=dippks(2);                                       % W2 intensity
        proc_var{record_no,86}=dipt(2);                                         % W2 time
        proc_var{record_no,87}=diparea(2);                                      % W2 area
        proc_var{record_no,88}=wri;                                             % WRI
        proc_var{record_no,89}=rhoc;                                            % rhoc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % Error trap for low SNR and message
        proc_var{record_no,1}=filename;                                         % filename
        proc_var{record_no,2}=metadata.datestring;                              % Date as text string
        proc_var{record_no,3}=metadata.bppvers;                                 % BP+ software version
        proc_var{record_no,4}=metadata.bppalgo;                                 % BP+ algorithm
        proc_var{record_no,5}=bRes_version;                                     % bRes version
        proc_var{record_no,6}=kres_v;                                           % kreservoir version
        proc_var{record_no,17}=metadata.snr;                                    % Signal to Noise Ratio
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
    xlsfile=strcat(folder_name, 'results\resdata.xls');
    % header = {'re_file' 're_basbp' 're_tbasbp' 're_minp' 're_intaopr' 're_maxaopr'...
    %     're_tmaxaopr' 're_intaoprlessdbp' 'date' 're_sam_rate'...
    %     're_intaoxsp' 're_maxaoxsp' 're_tmaxaoxsp' 're_aotn' 're_aopinf' 're_aopn'...
    %     're_aofita' 're_aofitb' 're_aorsq' 're_prob' 're_kres'...
    %     're_aitype' 're_hr' 're_sbp2' 're_intbapr' 're_maxbapr' 're_tmaxbapr'...
    %     're_intbaxsp' 're_maxbaxsp' 're_tmaxbap' 're_bafita' 're_bafitb'...
    %     're_barsq' 're_bapinf' 're_bapn' 're_ao_dpdt' 're_ba_dpdt'...
    %     're_pb_pf' 're_ri' 're_wf1i'  're_wf1t' 're_wf1a' 're_wbi' ...
    %     're_wbt' 're_wba' 're_wf2i'  're_wf2t' 're_wf2a'  're_wri' 're_rhoc' ...
    %     're_aosevr' 're_version' 're_quality' 'ao_ai'}; % header
    header = {'re_file' 're_date' 're_bppvers' 're_bppalgo' ...
        're_resvers' 're_kres' 're_sam_rate' 're_basbp' 're_tbasbp',...
        're_ba_dbp' 're_hr' 're_ba_map' 're_bapp' 're_aosbp' ...
        're_aodbp' 're_aopp' 're_snr' 're_rmssd' 're_sAI' 're_ppv'...
        're_rwttf' 're_rwttp' 're_sep' 're_quality' 're_ba_t1' 're_ba_p1'...
        're_ba_t2' 're_ba_p2' 're_pai' 're_ba_esp' 're_ba_dpdt' 're_ao_ed'...
        're_ao_esp' 're_ao_p1' 're_ao_p2' 're_aotr' 're_aitype' 're_ppar' ...
        're_ao_ap' 're_pmsys' 're_pmdia' 're_ao_tti' 're_ao_dti' 're_aosevr'...
        're_aodd' 're_ao_ai' 're_ai75' 're_aoti' 're_ao_ew' 're_ao_dpdt' ...
        're_pb' 're_pb_t' 're_pf' 're_pf_t' 're_PbPf' 're_ri' 're_intaopr'...
        're_maxaopr' 're_tmaxaopr' 're_intaoprlessdbp' 're_intaoxsp' ...
        're_maxaoxsp' 're_tmaxaoxsp' 're_aopinf' 're_aofita' 're_aofitb' ... 
        're_aorsq' 're_prob' 're_intbapr' 're_maxbapr' 're_tmaxbapr' ...
        're_intbaxsp' 're_maxbaxsp' 're_tmaxbap' 're_bafita' 're_bafitb' ...
        're_barsq' 're_bapinf' 're_wf1i' 're_wf1t' 're_wf1a' 're_wbi' ...
        're_wbt' 're_wba' 're_wf2i' 're_wf2t' 're_wf2a' 're_wri' ...
        're_rhoc'}; % header
   
    
    % % writetable
    Results_table=cell2table(proc_var, 'VariableNames',header);
    writetable(Results_table, xlsfile);
end
%% Tidy up
% rather than clearing the workspace for now I've left all non-redundant
% variables so that they can be used for debugging

% % message at end
msg=string(no_of_files);
msg=strcat(msg, " file(s) processed");
% msgbox(msg, 'Done');
disp(msg)

