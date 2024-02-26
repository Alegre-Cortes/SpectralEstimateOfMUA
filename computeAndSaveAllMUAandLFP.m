setParamsAndOptions;

%%
BaseDir = 'DataFolder';
RecordingDir = [BaseDir '\Record Node 101\experiment0\recording1'];
MessageDir = [RecordingDir '\events\Message_Center-904.0\TEXT_group_1'];

% This opens as 'memmapfile' the recorded raw signal.
% The script library 'analysis-tools-master' is required.
D = load_open_ephys_binary([RecordingDir '\structure.oebin'],'continuous',1,'mmap');

% Some required information to measure properly time and potentials.
SamplingRate = D.Header.sample_rate;            % in Hz.
RawSignal2uV = D.Header.channels(1).bit_volts;  % in uV to be multiplied to the raw signal.
RawSignalOffset = 0.0;                          % offset in uV to be subtracted to the raw signal.
RawSignalUnits = D.Header.channels(1).units;    % 'uV' string.
NoC = D.Header.num_channels;
t0 = double(D.Timestamps(1))/SamplingRate;      % time in seconds of the first sample.

ChannelSet = 1:NoC;
% ChannelSet = 1:2;

%% Loop on Channels...
%
tic
nc = 0;
for Channel = ChannelSet
   nc = nc + 1;
   fprintf('Ch. %d [%d/%d]...\n',Channel,nc,numel(ChannelSet));
   
   %% Load raw signal to analyze: to be tuned for the specific database...
   RawSig.value = (double(D.Data.Data.mapped(Channel,:))*RawSignal2uV - RawSignalOffset)/1000; % in mV
   RawSig.time = double(D.Timestamps)'/SamplingRate;
   RawSig.dt = diff(RawSig.time(1:2));

   %%  Compute baseline for MUA computation...
   %
   MovingWindowSize = 1/Options.LogMUA.FreqBand(1);
   DetrendingOrder = 0;
   [BaselineF, BaselinePSD] = plotMedianPSDofLFP(RawSig, MovingWindowSize, DetrendingOrder);
   close(gcf);

%    ndx = find(BaselineF>=Options.LogMUA.FreqBand(1) & BaselineF<=Options.LogMUA.FreqBand(2));
   ndx = find(BaselineF>Options.LogMUA.FreqBand(1)-0.1 & BaselineF<Options.LogMUA.FreqBand(2)+0.1);
%    ndxF = find(F>MUAFreqBand(1)-0.1 & F<MUAFreqBand(2)+0.1);
   BaselineF = BaselineF(ndx);
   BaselinePSD = BaselinePSD(ndx);
   
   %% Compute spectral estimate of MUA...
   %
   [MUA.time, MUA.value] = computeSpectralEstimateOfMUA(RawSig.time, RawSig.value, ...
      Options.LogMUA.FreqBand, 0, BaselinePSD);
   MUA.dt = diff(MUA.time(1:2));
   
   %% Compute smoothed log(MUA)...
   %
   logMUA.value = log(MUA.value);
   LowerBound = min(logMUA.value(logMUA.value>-Inf)); % Rectify null or negative values.
   logMUA.value(logMUA.value<=-Inf) = LowerBound;
   logMUA.time = MUA.time;
   logMUA.dt = MUA.dt;
   
   [logMUA.time, logMUA.value] = computeMovingAverageOfWF(logMUA, ...
      round(Options.LogMUA.SmoothingWindow / logMUA.dt), 1);
   logMUA.dt = diff(logMUA.time(1:2));
   
   %% Look for a reference value for logMUA.
   %  To be the final one, it should correspond to no firing activity...
   %
   FIRST_PEAK_METHOD = 1; % 'FirstPeak'
   LEFT_TAIL_METHOD = 2;  % 'LeftTail'
   
   MUARefMethod = FIRST_PEAK_METHOD;
   if isfield(Options.MUARef,'Method')
      if strcmpi(Options.MUARef.Method,'LeftTail')
         MUARefMethod = LEFT_TAIL_METHOD;
      end
   end
   
   switch(MUARefMethod)
      case LEFT_TAIL_METHOD
         BIN_NUM = 100;
         ValRange = [min(logMUA.value) max(logMUA.value)];
         X = linspace(ValRange(1), ValRange(2), BIN_NUM);
         Y = ksdensity(logMUA.value,X);
         idx = find(Y>=Options.MUARef.MinDensity,1,'first');
         logMUA.baseline = X(idx) + Options.MUARef.MUAShift;
         
      otherwise % FIRST_PEAK_METHOD
         % Computes the multimodal distribution of resampled MUA...
         ModeParams = plotMultimodalHistogram(logMUA.value);
         close(gcf);
         
         % Set MUA reference to the peak at smallest MUA...
         logMUA.baseline = ModeParams.Mu(1);
   end
   
   %% Detrend the computed logMUA baseline, and plot the histogram
   %  of log(MUA)...
   %
   logMUA.value = logMUA.value - logMUA.baseline;
   MUA.value = MUA.value * exp(logMUA.baseline);
   
   %% Store the result in the MEAMUA structure.
   %
%    MEAMUA.logsmoothed = true;
   MEAMUA.logsmoothed = false;
   if nc == 1
      if MEAMUA.logsmoothed
         MEAMUA.time = logMUA.time;
      else
         MEAMUA.time = MUA.time;
      end
      MEAMUA.dt = diff(MEAMUA.time(1:2));
      MEAMUA.values = nan(numel(ChannelSet),numel(MEAMUA.time));
      MEAMUA.channels = ChannelSet;
   end
   if MEAMUA.logsmoothed
      MEAMUA.values(nc,:) = logMUA.value;
   else
      MEAMUA.values(nc,:) = MUA.value;
   end
   
   %% Compute LFP, resampling and filtering the raw signal...
   %
   % Subsample LFP...
   WindowSize = round(1/RawSig.dt/Options.LFP.SubsamplingRate);
   [LFP.time,LFP.value] = computeMovingAverageOfWF(RawSig, WindowSize, WindowSize);
   LFP.dt = diff(LFP.time(1:2));
   
   % High-pass filtering of raw signal (LFP), if required...
   if isfield(Options.LFP,'HPFCutOffFreq')
      %             disp('LFP high-pass filtering...')
      [lpLFP.time,lpLFP.value] = computeMovingAverageOfWF(LFP, ...
         round(1/Options.LFP.HPFCutOffFreq/LFP.dt), ...
         round(1/Options.LFP.HPFCutOffFreq/LFP.dt/10));
      LFP.value = LFP.value - interp1(lpLFP.time,lpLFP.value,LFP.time,'pchip');
      clear lpLFP
   end
   
   %% Low-pass filtering of LFP...
   %
   [ssLFP.time,ssLFP.value] = computeMovingAverageOfWF(LFP, ...
      round(Options.LFP.SmoothingWindow/LFP.dt), 1);
   
   %% Store the result in the MEALFP structure.
   %
%    MEALFP.smoothed = true;
   MEALFP.smoothed = false;
   if nc == 1
      if MEALFP.smoothed
         MEALFP.time = ssLFP.time;
      else
         MEALFP.time = LFP.time;
      end
      MEALFP.dt = diff(MEALFP.time(1:2));
      MEALFP.values = nan(numel(ChannelSet),numel(MEALFP.time));
      MEALFP.channels = ChannelSet;
   end
   if MEALFP.smoothed
      MEALFP.values(nc,:) = ssLFP.value;
   else
      MEALFP.values(nc,:) = LFP.value;
   end

end
toc

save('MEAMUALFP.mat','MEAMUA','MEALFP','Options','-v7.3');
