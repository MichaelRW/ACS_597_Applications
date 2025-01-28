%% Plot effective NR for Hearing Protector vs time
% Carl Howard 22/9/2016
%
% Calculate the effective noise reduction as a percentage of the time the
% hearing protection device is worn, for HPD that provide 10, 15, 20, 25,
% 30 dB of noise reduction.

% The hearing protector's rated Noise Reduction
NR_HPD = [10 15 20 25 30];  % dB

% The ambient noise level in the environment
SPL1 = 100;     % dB - the exposure level

% Assume an 8 hour day, 85 dBA over 8 hours.
minutes_8hr = 8*60;     % Number of minutes in 8 hours

% Number of minutes hearing protection is NOT worn
minutes_HPD_off = [0:1:120];    % minutes

% number of minutes hearing protection IS worn
minutes_HPD_on = minutes_8hr-minutes_HPD_off;

% Do a for loop to generate the curves for each of the HPD ratings
for jj=1:length(NR_HPD),
    Leq2(jj,:) = 10*log10( (minutes_HPD_off/minutes_8hr)*10^((SPL1)/10) +  ...
                 (minutes_HPD_on/minutes_8hr)*10^((SPL1-NR_HPD(jj))/10) );
    NR_effective(jj,:) = SPL1 - Leq2(jj,:);
end;

% Plot the graph
p1=plot(minutes_HPD_on/minutes_8hr*100,NR_effective,'k-','linewidth',2.0);
set(gca,'fontsize',14);
set(gca, 'xdir','reverse')
title('Effective Noise Reduction of HPD','FontWeight','Normal');
xlabel({'Percentage of time HPD Worn';'During Noise Exposure [%]'});
ylabel('Effective NR of HPD [dB]');

grid on
hA=gca;
set(hA,'GridAlpha',0.8);       % make gridlines more visible
set(hA,'XMinorTick','on','YMinorTick','on')
hA.YAxis.MinorTickValues = [0:1:40];
hA.XAxis.MinorTickValues = [75:1:100];


%% Example 1
% person has a HPD with nominal 30 dB Noise reduction
% but only wears for 95% of the time. What is the effective NR?
t_off = 24;     % min time not worn
t_shift = 8*60;   % min time of noise exposure
percent_worn = ((t_shift-t_off)/t_shift) * 100
ex1_Lp2=10*log10( ((t_shift-t_off)/t_shift)*10^((SPL1-30)/10)+    ...
        (t_off/t_shift)*10^(SPL1/10));
ex1_NR_eff=SPL1-ex1_Lp2

hold on
p2=plot([percent_worn percent_worn 100 ],   ...
        [0 ex1_NR_eff ex1_NR_eff],'k:','linewidth',2);