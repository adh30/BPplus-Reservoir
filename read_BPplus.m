function  [data, metadata, ss, ba, ao]=read_BPplus(folder_name, filename, Npoly, Frame)
% version 1.1 [ADH]
% fixed bugs or vulnerabilities in earlier version
% - Out-of-bounds indexing at start of ba.p_all
% - Out-of-bounds indexing at end of ba.p_all
% - Out-of-bounds indexing at start of ao.p_all
% - Out-of-bounds indexing at end of ao.p_all
% - Out-of-bounds indexing when calculating pulse waveforms
% - Insufficient pulse indices before diff() and allocation
% STILL NEED TO CHECK WHETHER CALIBRATION IS CORRECT - assumes that ss.p_av is zero referenced and in mmHg. Compare with stated SBP and DBP from conventional measurement.
%%
if ~isfile([folder_name filename])
    disp("No BP+ xml file exists")
    return;
end
[data] = xml2struct([folder_name filename]);

% Read Data from either CardioScope legacy xml or BPplus xml
if sum(strcmp(fieldnames(data), 'CardioScope')) == 1
    [metadata, ba, ao, ss] = read_BPplusCardioScope(data, Npoly, Frame);
else
    [metadata, ba, ao, ss] = read_BPplusBPplus(data);
end

%% categorise quality based on SNR
if metadata.snr>=12
    metadata.quality='Excellent';
elseif metadata.snr>=9
    metadata.quality='Good';
elseif metadata.snr>=6
    metadata.quality='Acceptable';
elseif metadata.snr>0
    metadata.quality='Poor';
else
    metadata.quality='Unacceptable';
end

% do we extract the measurement ID number?
%file ID (i.e. filname without the .xml
metadata.fileID = extractBefore( filename , '.xml');

%% brachial pulses
% replace values at start and end <ba.dbp and >SBP with NaN
% deal with low early values
for i = 1:min(200,length(ba.p_all))
    if ba.p_all(i)<ba.dbp
        ba.p_all(i) = NaN;
    end
end
% deal with high end values
for i = max(1,length(ba.p_all)-200):length(ba.p_all)
    if ba.p_all(i)>ba.sbp
        ba.p_all(i) = NaN;
    end
end

% remove NaNs (If removed, must adjust offsets)
%ba.p_all = ba.p_all(~isnan(ba.p_all)); this will break offsets to start of pulses ...

% calculate individual selected pulses for later display
numgoodpulses = length(ss.SelectedPulseIndexes)-1;
pulseindex = int32(ss.pulseStartIndexes);

if numel(pulseindex) < 2
    error('read_BPplus:InsufficientPulseIndexes', ...
          'At least two pulse start indexes are required');
end

pulselengths = diff(pulseindex);
ss.pulsewaveforms = NaN(max(pulselengths),numgoodpulses);

for i = 1:numgoodpulses
    crop=ba.p_all(pulseindex(i):pulseindex(i+1));
    ss.pulsewaveforms(1:length(crop),i)=crop;
end

% TODO: brachial average beat or supra-systolic average beat?
ba.p_av=ss.p_av;
% error trap if denominator is zero
denom = max(ba.p_av)-min(ba.p_av);

if denom == 0
    error('read_BPplus:ZeroPulseAmplitude', ...
          'Average pulse waveform has zero amplitude');
end

ss_cal_p = ba.pp/denom;
ba.p_av=ba.dbp+(ba.p_av*ss_cal_p);                                             
ss.dpdt=ss.dpdt*ss_cal_p;                                                      % correcting to mmHg/s

% aortic pulses (If removed, must adjust offsets)
%ao.p_all=ao.p_all(ao.pulseStartIndexes(1):ao.pulseStartIndexes(end));

% replace values at start and end <ba.dbp and >SBP with NaN
% deal with low early values
for i = 1:min(200,length(ao.p_all))
    if ao.p_all(i)<ba.dbp
        ao.p_all(i) = NaN;
    end
end
% deal with high end values
for i = max(1,length(ao.p_all)-200):length(ao.p_all)
    if ao.p_all(i)>ao.sbp
        ao.p_all(i) = NaN;
    end
end
% % remove NaNs (If removed, must adjust offsets)
% ao.p_all = ao.p_all(~isnan(ao.p_all));

end
