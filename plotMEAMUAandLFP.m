%% Set white-orange-purple colormap for MUA and blue-red colormap for LFP.
Purple = [hex2dec('95')/255 hex2dec('63')/255 hex2dec('94')/255];
Orange = [hex2dec('f4')/255 hex2dec('95')/255 hex2dec('05')/255];
LightOrange = (Orange+1)/2;
White = [1 1 1];
% CM = multigradientColormap([White; Orange; Purple]);
MUACM = multigradientColormap([White; LightOrange; Orange; Purple]);
% CM = multigradientColormap([White; LightOrange; Purple]);

LFPCM = gradedColormap([0 0 1],[1 0 0]);


%% MUA rasterplot
MUARange = [-0.5 2.0];

figure

ax1 = subplot(4,1,1);
plot(MEAMUA.time,mean(MEAMUA.values),'k','LineWidth',0.75)
xlim(MEAMUA.time([1 end]))
ylim(MUARange)
hcb = colorbar();
hcb.Visible = 'off';
set(gca,'Layer','top','TickDir','out','Box','on','YDir','normal','XTickLabels',{})
ylabel('MEA log(MUA)')
   
ax2 = subplot(4,1,2:4);
imagesc(MEAMUA.time,1:size(MEAMUA.values),MEAMUA.values,MUARange)
colormap(MUACM);
hcb = colorbar();
hcb.Label.String = 'log(MUA)';
set(gca,'Layer','top','TickDir','out','Box','on','YDir','normal')
xlim(MEAMUA.time([1 end]))
xlabel('Time [s]')
ylabel('Channels')

linkaxes([ax1 ax2],'x')

FigSize = [5 4];
FigName = 'MEAMUAvsTime';
set(gcf,'PaperUnit','inch','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf',sprintf('%s.pdf',FigName),'-painters')
% print('-dpng',sprintf('%s.png',FigName),'-r300')

%% LFP rasterplot
LFPRange = [-1 +1]*0.3;

figure

ax1 = subplot(4,1,1);
plot(MEALFP.time,mean(MEALFP.values),'k','LineWidth',0.75)
xlim(MEALFP.time([1 end]))
ylim(LFPRange)
hcb = colorbar();
hcb.Visible = 'off';
set(gca,'Layer','top','TickDir','out','Box','on','YDir','normal','XTickLabels',{})
ylabel('MEA LFP [mV]')
   
ax2 = subplot(4,1,2:4);
imagesc(MEALFP.time,1:size(MEALFP.values),MEALFP.values,LFPRange)
colormap(LFPCM);
hcb = colorbar();
hcb.Label.String = 'LFP [mV]';
set(gca,'Layer','top','TickDir','out','Box','on','YDir','normal')
xlim(MEALFP.time([1 end]))
xlabel('Time [s]')
ylabel('Channels')

linkaxes([ax1 ax2],'x')

FigSize = [5 4];
FigName = 'MEALFPvsTime';
set(gcf,'PaperUnit','inch','PaperPosition',[0 0 FigSize],'PaperSize',FigSize);
print('-dpdf',sprintf('%s.pdf',FigName),'-painters')
% print('-dpng',sprintf('%s.png',FigName),'-r300')
