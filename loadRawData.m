% This is an example code to extract the raw signal collected from the two
% MEA implanted in Bart, each containing 96 recording channels.

BaseDir = 'D:\Home\Mattia\Qsync\Research\ExperimentalData\2022.AnaestisiaFadeOut_Bart.Ferraina';
RecordingDir = [BaseDir '\Record Node 101\experiment0\recording1'];
MessageDir = [RecordingDir '\events\Message_Center-904.0\TEXT_group_1'];

% This opens as 'memmapfile' the recorded raw signal.
% The script library 'analysis-tools-master' is required.
D = load_open_ephys_binary([RecordingDir '\structure.oebin'],'continuous',1,'mmap')

% Some required information to measure properly time and potentials.
SamplingRate = D.Header.sample_rate;            % in Hz.
RawSignal2uV = D.Header.channels(1).bit_volts;  % in uV to be multiplied to the raw signal.
RawSignalOffset = 0.0;                          % offset in uV to be subtracted to the raw signal.
RawSignalUnits = D.Header.channels(1).units;    % 'uV' string.
NoC = D.Header.num_channels;
t0 = double(D.Timestamps(1))/SamplingRate;      % time in seconds of the first sample.

%% Extract a data chunk and plot it.
% TimeRange = 9200 + [0 10]; % This a time window of 2 seconds following anaesthesia injection.
% ndx = floor((TimeRange - t0)*SamplingRate + 1);
% ndx = ndx(1):ndx(end);

ndx = 1:(SamplingRate/60):numel(D.Timestamps);

LFP = (double(D.Data.Data.mapped(:,ndx))*RawSignal2uV - RawSignalOffset)/1000; % in mV
% LFPrange = max(abs(LFP(:)))*[-1 1]/2;
LFPrange = 0.3*[-1 1];

figure
imagesc((ndx-ndx(1))/SamplingRate,1:NoC,LFP,LFPrange)
colorbar
xlabel('Time [s]')
ylabel('Channels')
title(sprintf('LFP [mV], t_0 = %.4g',(ndx(1)-1)/SamplingRate+t0))