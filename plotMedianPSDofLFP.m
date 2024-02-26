function [F, PSDMedian] = plotMedianPSDofLFP(LFP, MovingWindowSize, DetrendingOrder)
%
%  [F, PSDMedian] = plotMedianPSDofLFP(LFP, MovingWindowSize[, DetrendingOrder])
%
%   Copyright 2014 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Jan. 10, 2014
%

if exist('DetrendingOrder','var') == 0
   DetrendingOrder = 2; % Remove linear trend...
end

%---
%  Arranges WFs for spectral analysis...
%---
SamplingFreq = 1/diff(LFP.time(1:2));
WindowSize = round(MovingWindowSize*SamplingFreq);
SampleNum = floor(length(LFP.value)/WindowSize);

LFPs = reshape(LFP.value(1:SampleNum*WindowSize), WindowSize, SampleNum);

%---
% Polinomial detrending...
%---
if DetrendingOrder > 0 % 0-th order detrending...
   LFPs = LFPs - repmat(mean(LFPs), WindowSize, 1); 
end
if DetrendingOrder > 1 % 1-st order detrending...
   Xx = repmat(linspace(-WindowSize/2,WindowSize/2,WindowSize)', 1, sn);
   LFPs = LFPs -  Xx .* repmat(mean(diff(LFPs,1,1)), WindowSize, 1);
end
if DetrendingOrder > 2 % 2-nd order detrending...
   LFPs = LFPs - 1/2 * (Xx.^2) .* repmat(mean(diff(LFPs,2,1)), WindowSize, 1);
end
if DetrendingOrder > 3 % 3-rd order detrending...
   LFPs = LFPs - 1/6 * (Xx.^3) .* repmat(mean(diff(LFPs,3,1)), WindowSize, 1);
end


%---
%  Spectral estimates of WFs...
%---
Y = fft(LFPs, WindowSize);
Pyy =  Y.* conj(Y) / WindowSize;
% Pyy =  Y.* conj(Y) / WindowSize / SamplingFreq;
F = SamplingFreq * (0:WindowSize/2)' / WindowSize;

ndxF = 2:length(F);
F = F(ndxF);
Pyy = Pyy(ndxF,:);

PSDMedian = prctile(Pyy', 50)';

%---
%   Plot average power spectrum and its distribution...
%---
figure
hold on

% Patches representing different percetiles...
XP = [F; flipud(F)]';
for prc = 10:10:40
   clr = [1-prc/100 1-prc/100 1];
   YP = prctile(Pyy', [prc 100-prc]);
   patch(XP, [YP(1,:) fliplr(YP(2,:))], clr, 'EdgeColor', 'none');
end

plot(F, PSDMedian, '.-b', 'LineWidth', 1.);

set(gca, 'XLim', [F(1) F(end)], 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca, 'TickDir', 'out', 'Layer', 'top', 'Box', 'off');

xlabel('\omega/2\pi (Hz)');
ylabel('LFP p.s.d. (\muV^2/Hz)');

% set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
% print('-deps2c', 'MedianPSDofLFP.eps');
