function ModeParams = plotMultimodalHistogram(Values, ValRange, BinNum)
%
%  ModeParams = plotMultimodalHistogram(Values[, ValRange[, BinNum]])
%
%  Makes use of the Statistics toolbox. The histogram is smoothed with
%  ksdensity(). The fit is performed looking at a suited mixture of
%  Gaussian distributions.
%  
%  Maurizio Mattia @ 2013, ver 1.0
%

GRAY_CLR = repmat(0.9, 1, 3);
BIN_NUM = 100;
MIN_RESIDUAL_STD = 0.05; % 5% of the maximum density estimated.


%% Parameter settings...
%
if exist('BinNum','var')
   if BinNum>0
      BIN_NUM = BinNum;
   end
end
if exist('ValRange','var') ~= 1
   ValRange = [min(Values) max(Values)];
else
   if isempty(ValRange)
      ValRange = [min(Values) max(Values)];
   end
end
if size(Values,1) < size(Values,2)
   Values = Values';
end
X = linspace(ValRange(1), ValRange(2), BIN_NUM);


%% Fit different Gaussian mixtures...
%
Y = ksdensity(Values,X);

MeanValue = mean(Values);
ndxX = find(X <= MeanValue);

Options = statset('MaxIter',1000);
% FitIsBad = 1;
% k = 1;
% while FitIsBad && k<=4
%    obj = gmdistribution.fit(Values,k,'Options',Options);
%    objY = pdf(obj,X')';
%    if std(objY(ndxX)-Y(ndxX))/max(Y) < MIN_RESIDUAL_STD
%       FitIsBad = 0;
%    else
%       k = k + 1;
%    end
% end
MaxNoG = 6;
AIC = zeros(1,MaxNoG);
for k = 1:MaxNoG
   obj = gmdistribution.fit(Values,k,'Options',Options);
   AIC(k) = obj.AIC;
end
DeltaAIC = abs(diff(AIC));
% k = find(DeltaAIC<0.005,1,'first');
[~,k] = min(DeltaAIC);
obj = gmdistribution.fit(Values,k,'Options',Options);

%% Compose fit results for the output...
%
ModeParams.Mu = obj.mu;
ModeParams.Sigma = squeeze(obj.Sigma);
ModeParams.Ampl = obj.PComponents';
[ModeParams.Mu,ndx] = sort(ModeParams.Mu);
ModeParams.Sigma = ModeParams.Sigma(ndx);
ModeParams.Ampl = ModeParams.Ampl(ndx);


%% Plot resulting fit...
%
figure
hold on

YP = [Y zeros(size(Y))];
XP = [X fliplr(X)];
patch(XP,YP,GRAY_CLR,'EdgeColor','none');
plot(X,Y,'k','LineWidth',1.)

% THIS IS WRONG BECAUSE GAUSSIAN MIXTURE IS MULTIPLICATIVE, NOT ADDITIVE...
% clr = copper(obj.NComponents+1);
% for k = 1:obj.NComponents
%    plot(X,mygauss([ModeParams.Mu(k) ModeParams.Sigma(k) ModeParams.Ampl(k)],X)/(ModeParams.Sigma(k)*sqrt(2*pi)),'-','Color',clr(k+1,:));
% end

objY = pdf(obj,X')';
plot(X,objY,'g-','LineWidth',1.);
YRange = get(gca,'YLim');

clr = copper(obj.NComponents+2);
for k = 1:obj.NComponents
   plot([0 0]+ModeParams.Mu(k),YRange,'--','Color',clr(k+1,:));
   ndx = find(X>=ModeParams.Mu(k),1,'first');
   text(X(ndx),max([objY(ndx) Y(ndx)])+diff(YRange)/20,[' ' num2str(ModeParams.Mu(k),3)],'Color',clr(k+1,:));
end
set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');
set(gca, 'XLim', ValRange);

xlabel('Value (a.u)');
ylabel('Density');
