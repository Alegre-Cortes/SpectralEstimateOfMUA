Options.PeriodToAnalyze = [0 Inf];

Options.SaveMUA = 0;

Options.LogMUA.FreqBand = [200 1500];
Options.LogMUA.SmoothingWindow = 0.040; % It depends on the sampling rate.

Options.LFP.SmoothingWindow = Options.LogMUA.SmoothingWindow; % s...
Options.LFP.SubsamplingRate = Options.LogMUA.FreqBand(1);   % Hz...
Options.LFP.HPFCutOffFreq = 0.2; % Hz...
Options.LFP.xWavScale = 10^6; % Conversion factor from V to uV (if xWav are in V).

Options.UpTimeLags.MaxAbsTimeLag = 0.250;         % Maximum reasonable time lag between electrodes...

Options.LFPMUAplot.MUArange = [-0.75 3.5];
Options.LFPMUAplot.LFPrange = [-400 400]; % for LFP

Options.MUARef.Method = 'LeftTail'; % Or 'FirstPeak'
Options.MUARef.MinDensity = 0.05;
Options.MUARef.MUAShift = 0.35;

PrintDev = '-dpdf';

%% Useful colors and palettes.
%
Red = [1 0 0];
Green = [0 1 0];
Blue = [0 0 1];
Yellow = Red + Green;
YellowOrange = 0.75*Yellow + 0.25*Red;
Orange = 0.5*Yellow + 0.5*Red;
BluePurple = 0.75*Blue + 0.25*Red;

RedBlueCB = gradedColormap(Blue,Red);
