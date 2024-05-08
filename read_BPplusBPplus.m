function [metadata, ba, ao, ss] = read_BPplusBPplus(data)
    %% Read Data from BPplus xml
    result = data.BPplus.Results.Result;
    measDataLogger = data.BPplus.MeasDataLogger;
    % meta data
    metadata.bppvers=measDataLogger.Attributes.version;                                 % Software version
    if isfield(measDataLogger.Attributes,'nibp')
        metadata.nibp=measDataLogger.Attributes.nibp;                                   % nibp type
    end
    metadata.datestring=measDataLogger.Attributes.datetime;                             % Date as text string
    metadata.guid = measDataLogger.Attributes.guid;                                     % Unique ID for this measurement
    metadata.bppalgo=result.Attributes.algorithm_revision;                              % Software algorithm
    metadata.samplerate=str2double(measDataLogger.SampleRate.Text);                     % sample rate, Hz
    metadata.snr=str2double(result.SNR.Text);                                           % Signal to noise ratio, dB
    metadata.patienet_id = '';
    metadata.notes = '';

    metadata.RawSuprasystolicPressure=measDataLogger.RawSuprasystolicPressure.Text;     % raw base 64 data

    % Mode is only available in later versions of xml
    if isfield(data,'NibpMode')
        metadata.mode = string(data.BPplus.MeasDataLogger.NibpModeUsed.Text);           % Measurement mode
    end

    ss.prv=str2double(result.sPRV.Text);                                                % RMSSD of beats in suprasystolic signal
    ss.ai=str2double(result.sAI.Text);                                                  % AI from suprasystolic signal
    ss.dpdt=str2double(result.sDpDtMax.Text);                                           % dp/dt from suprasystolic signal in uncorrected units
    %ss.harm= Can be calculated from sAI if needed
    ss.pp=str2double(result.sPP.Text);                                                  % Suprasystolic Pulse Pressure in the Cuff, mmHg (PP in cuff is small)
    ss.ppv=str2double(result.sPPV.Text);                                                % Pulse pressure variation, % [?]
    ss.RWTTFoot=str2double(result.sRWTTFoot.Text);                                      % reflected wave transit time from foot of suprasystolic signal
    ss.RWTTPeak=str2double(result.sRWTTPeak.Text);                                      % reflected wave transit time from peak of suprasystolic signal
    ss.sep=str2double(result.sSEP.Text)/1000;                                           % Systolic ejection period, s
    ss.Tn=split(result.sAveragePulsePointsIndexes.Text,',');                            % Times of pulse points (1-6) see diagram. time=Tn/SampleRate

    %% suprasystolic rhythm, average beat & start of pulses.
    ss.p_all=str2double(split(result.sBaseLined.Text,','));                             % not scaled to arterial mmHg
    ss.pulseStartIndexes=str2double(split(result.sPulseStartIndexes.Text,','))+1;
    ss.SelectedPulseIndexes=str2double(split(result.sSelectedPulseIndexes.Text,','))+1; % +1 as MATLAB is not zero indexed.
    ss.p_av = str2double(split(result.sAveragePulse.Text,','));                         % not scaled to arterial mmHg
    ss.averagePulsePointsIndexes=str2double(split(result.sAveragePulsePointsIndexes.Text,','))+1;

    % nibp brachial measurement
    % FIXME: naming confusion between NIBP values and those from the calculated ba waveforms.
    ba.sbp=str2double(measDataLogger.Sys.Text);                                         % brachial systolic BP (NIBP), mmHg
    ba.dbp=str2double(measDataLogger.Dia.Text);                                         % brachial diastolic BP (NIBP), mmHg
    ba.pp=ba.sbp-ba.dbp;                                                                % brachial pulse pressure (NIBP), mmHg
    ba.map=str2double(measDataLogger.Map.Text);                                         % brachial mean arterial pressure (NIBP), mmHg
    ba.hr=str2double(measDataLogger.Pr.Text);                                           % heart rate (NIBP), bpm
    %% estimated brachial rhythm
    ba.p_all=str2double(split(result.baEstimate.Text,','));

    % calculated aortic measurements
    ao.sbp=str2double(result.cSys.Text);                                                % cSBP calculated by BP+, mmHg
    ao.dbp=str2double(result.cDia.Text);                                                % cDBP calculated by BP+, mmHg
    ao.map=str2double(result.cMap.Text);                                                % cMAP calculated by BP+, mmHg
    ao.pp=ao.sbp-ao.dbp;                                                                % cPP calculated by BP+, mmHg
    if isfield(result,'cST')
        ao.ed=str2double(result.cST.Text)/1000;                                         % cED calculated by BP+. duration of systole. NOTE: BP+ cED is a %
    else
        ao.ed=[];
    end
    
    %% aortic rhythm, average beat & start of pulses.
    ao.p_all=str2double(split(result.cEstimate.Text,','));

    % assumed delay from aortic pulse to cuff.
    ao.lag=str2double(result.cLag.Text);
    ao.pulseStartIndexes=str2double(split(result.cPulseStartIndexes.Text,','))+1;
    ao.p_av=str2double(split(result.cAveragePulse.Text,','))';
    if isfield(result,'cAveragePulsePointsIndexes')
        ao.averagePulsePointsIndexes=str2double(split(result.cAveragePulsePointsIndexes.Text,','))+1;
    end
    if isfield(result,'cAveragePulsePointsIndexes')
        ao.averagePulsePointsIndexes=str2double(split(result.cPulseStartIndexes.Text,','))+1;
    end

    %pPX brachial pulsatility index (PP/ba.map)
    %cPX aortic pulsatility index (aoPP/ba.map)
end