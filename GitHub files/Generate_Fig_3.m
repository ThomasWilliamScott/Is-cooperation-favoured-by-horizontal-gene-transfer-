% Script to generate figure 3 (plot of effective plasmid transfer 
% frequency).

close all
clearvars
clc

hold on % allows us to plot multiple lines on the same graph

plasfreq=[0:0.001:1]; % defines the range of plasmid frequency values to be used on the x axis.
n=4; % number of founder cells
beta=0.95; % plasmid transfer probability

% evaluate expression for effective plasmid transfer frequency.
realised = (1-plasfreq) .* ((n-1)/n) .* beta; 

% plot reference line (effective plasmid transfer frequency) in black
plot(plasfreq,realised,'LineWidth',2,'Color','k')

% change plasmid transfer probability for second line
beta=0.6;

% re-evaluate expression for effective plasmid transfer frequency for second line.
realised = (1-plasfreq) .* ((n-1)/n) .* beta;

% plot second line (reduced plasmid transfer probability)
plot(plasfreq,realised,'LineWidth',2,'LineStyle',':','Color',[0.8500, 0.3250, 0.0980]	)

% change plasmid transfer probability back to reference value
beta=0.95;

% change number of founder cells for third line
n=2;

% re-evaluate expression for effective plasmid transfer frequency for third line.
realised = (1-plasfreq) .* ((n-1)/n) .* beta;

% plot third line (reduced founder cells)
plot(plasfreq,realised,'LineWidth',2,'LineStyle',':','Color',[0.9290, 0.6940, 0.1250]	)

% make y axis vary between zero and one.
ylim([0 1])

% formatting
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
yticks([0:0.2:1])
