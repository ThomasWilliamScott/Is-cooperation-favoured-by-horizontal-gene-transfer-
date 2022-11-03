% Script to generate supplementary figure 1 (plot of relatedness function). 

close all
clearvars
clc

hold on % allows us to plot multiple lines on the same graph

pcR=[0:0.001:1]; % defines the range of cooperative plasmid frequency values to be used on the x axis.
pd=0; % defector plasmid frequency value
n=4; % number of founder cells
beta=0.95; % plasmid transfer probability

% plasmid relatedness function
R =((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - (pcR + pcR.^2.*(n - 1) - (2.*beta.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n - (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2 + (beta.^2.*pcR.*(n - 1).*(n - 2).*(pcR + pd - 1).^2.*(n.*pcR - 3.*pcR + 1))./n.^2)./n)./((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - pcR + (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2);

% plot first line (reference line) in black
plot(pcR,R,'LineWidth',2,'Color','k')

% plot line showing chromosomal relatedness (given by 1/N)
yline(1/n,'--','LineWidth',2)

% change plasmid transfer probability for second line
beta=0.6;

% re-evaluate plasmid relatedness function (for second line)
R =((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - (pcR + pcR.^2.*(n - 1) - (2.*beta.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n - (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2 + (beta.^2.*pcR.*(n - 1).*(n - 2).*(pcR + pd - 1).^2.*(n.*pcR - 3.*pcR + 1))./n.^2)./n)./((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - pcR + (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2);

% plot second line (reduced plasmid transfer probability)
plot(pcR,R,'LineWidth',2,'LineStyle',':','Color',[0.8500, 0.3250, 0.0980]	)

% change plasmid transfer probability back to the reference value
beta=0.95;

% change number of founder cells for third line
n=2;

% re-evaluate plasmid relatedness function (for third line)
R =((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - (pcR + pcR.^2.*(n - 1) - (2.*beta.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n - (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2 + (beta.^2.*pcR.*(n - 1).*(n - 2).*(pcR + pd - 1).^2.*(n.*pcR - 3.*pcR + 1))./n.^2)./n)./((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - pcR + (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2);

% plot third line (reduced founder cells)
plot(pcR,R,'LineWidth',2,'LineStyle',':','Color',[0.9290, 0.6940, 0.1250]	)

% change number of founder cells back to reference value
n=4;

% change defector plasmid frequency for fourth line
pd=0.6; 

% re-evaluate plasmid relatedness function (for fourth line)
R =((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - (pcR + pcR.^2.*(n - 1) - (2.*beta.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n - (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2 + (beta.^2.*pcR.*(n - 1).*(n - 2).*(pcR + pd - 1).^2.*(n.*pcR - 3.*pcR + 1))./n.^2)./n)./((pcR - (beta.*pcR.*(n - 1).*(pcR + pd - 1))./n).^2 - pcR + (beta.^2.*pcR.*(n - 1).*(pcR + pd - 1).*(n.*pcR - 2.*pcR + 1))./n.^2);

% plot fourth line (increased defector plasmid frequency)
plot(pcR(pcR<=(1-pd)),R(pcR<=(1-pd)),'LineWidth',2,'LineStyle',':','Color',[0.4940, 0.1840, 0.5560]	)

% make y axis vary between zero and one.
ylim([0 1])

% formatting
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
yticks([0:0.2:1])