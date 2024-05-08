function [metadata, ba, ao, ss] = read_BPplusCardioScope(data, Npoly, Frame)
    % extract CardioScope values
    result = data.CardioScope.Results.Result;
    measDataLogger = data.CardioScope.MeasDataLogger;
    % meta data
    metadata.bppvers=measDataLogger.Attributes.version;                                 % Software version
    metadata.nibp=measDataLogger.Attributes.nibp;                                       % nibp type
    metadata.datestring=measDataLogger.Attributes.datetime;                             % Date as text string
    metadata.guid = measDataLogger.Attributes.device_id;                                     % Unique ID for this measurement
    metadata.bppalgo=result.Attributes.algorithm_revision;                              % Software algorithm
    metadata.samplerate=str2double(measDataLogger.SampleRate.Text);                     % sample rate, Hz
    metadata.snr=str2double(result.SNR.Text);                                           % Signal to noise ratio, dB
    metadata.patient_id = '';                                                           
    metadata.notes = '';

    metadata.RawSuprasystolicPressure=measDataLogger.SSBuff.Text;                       % raw base 64 data

    ss.prv=str2double(result.RMSSD.Text);                                               % RMSSD of beats in suprasystolic signal
    ss.ai=str2double(result.ssAI.Text);                                                 % AI from suprasystolic signal
    ss.dpdt=str2double(result.ssDpDtMax.Text);                                          % dp/dt from suprasystolic signal in uncorrected units
    %ss.harm=str2double(result.ssHARM.Text);                                            % Normalized sAI
    ss.pp=str2double(result.ssPP.Text);                                                 % Suprasystolic Pulse Pressure
    ss.ppv=str2double(result.ssPPV.Text);                                               % Suprasystolic Pulse Pressure Variation, %
    ss.RWTTFoot=str2double(result.ssRWTTFoot.Text);                                     % reflected wave transit time from foot of suprasystolic signal
    ss.RWTTPeak=str2double(result.ssRWTTPeak.Text);                                     % reflected wave transit time from peak of suprasystolic signal
    ss.sep=str2double(result.ssSEP.Text);                                               % Systolic ejection period, s
    ss.Tn=split(result.ssAverageBeatPointsIdxs.Text,',');                               % Offset to characteristic points (Tn/SampleRate is time)

    %% suprasystolic rhythm, average beat & start of pulses.
    ss.p_all=[];                                                                        % CardioScope did not save sBaselinddData. TODO calculate it?
    ss.pulseStartIndexes=str2double(split(result.ssBeatStartIdxs.Text,',')) + 1;
    ss.SelectedPulseIndexes=  (0:1:length(ss.pulseStartIndexes)-1)+1;                   % not saved in XML. auto select all for now?
    ss.p_av = str2double(split(result.ssAverageBeat.Text,','));                         % not scaled to arterial mmHg
    ss.averagePulsePointsIndexes=str2double(split(result.ssAverageBeatPointsIdxs.Text,','))+1;

    % nibp brachial measurement
    % FIXME: naming confusion between NIBP values and those from the calculated ba waveforms.
    ba.sbp=str2double(measDataLogger.Sys.Text);                                         % brachial systolic BP (NIBP), mmHg
    ba.dbp=str2double(measDataLogger.Dia.Text);                                         % brachial diastolic BP (NIBP), mmHg
    ba.pp=ba.sbp-ba.dbp;                                                                % brachial pulse pressure (NIBP), mmHg
    ba.map=str2double(measDataLogger.Mean.Text);                                        % mean arterial pressure (NIBP), mmHg
    ba.hr=str2double(measDataLogger.Hr.Text);                                           % heart rate (NIBP), bpm
    %% estimated brachial rhythm
    ba.p_all=str2double(split(result.baEstimate.Text,','));

    % calculated aortic measurements
    ao.sbp=str2double(result.aoSys.Text);                                               % cSBP calculated by BP+, mmHg
    ao.dbp=str2double(result.aoDia.Text);                                               % cDBP calculated by BP+, mmHg
    %ao.map=str2double(result.aoMap.Text);                                               % TODO: calculate from average beat? cMAP calculated by BP+, mmHg
    ao.pp=ao.sbp-ao.dbp;                                                                % cPP calculated by BP+, mmHg
    ao.ed = -1;                                                                          % TODO Calcualte ED as duration of systole. NOTE: BP+ cED is a %

    %% aortic rhythm, average beat & start of pulses.
    ao.p_all=str2double(split(result.aoEstimate.Text,','));

    % assumed delay from aortic pulse to cuff.
    ao.lag = 0.18;
    if (contains(data.CardioScope.MeasDataLogger.Attributes.software_version,"038"))
        ao.lag = 0.06;
    end

    % Calculate offset to start of central pulses from suprasystolic pulse start indexes and aoLag.
    ao.pulseStartIndexes = ss.pulseStartIndexes - round(ao.lag * metadata.samplerate);
    ao.pulseStartIndexes(ao.pulseStartIndexes < 1) = 1;

    ao.p_av = str2double(split(result.aoAverageBeat.Text,','));
    % not avilable ao.averagePulsePointsIndexes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIXME repartition and calc average pulse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% aortic average beat, ao_p_av
    a0=str2double(split(result.aoAverageBeat.Text,','));
    % create a double beat' to deal with the errors in definition of dbp
    a=[a0; a0];
    %plot(a);
    % identify start and end of beat as minima
    [~, locmin1]=min(a);
    [~, locmin2]=min(a(locmin1+50:end));
    locmin2=locmin2+locmin1+50;
    %plot(a(locmin1:locmin2));

    % filter derivative of new beat with SG to get rid of kinks at or around join
    aa = sgolayfilt(diff(a),Npoly,Frame); % filter the derivative - order and framelen defined above
    a1=cumsum(aa)-min(cumsum(aa))+min(a);    % reconstruct
    ao.p_av =a1(locmin1:locmin2-1)';         % crop to cycle
    clear a0 a aa a1 locmin1 locmin2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end