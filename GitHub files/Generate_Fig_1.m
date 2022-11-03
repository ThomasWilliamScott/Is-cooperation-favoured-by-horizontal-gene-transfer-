% Script to generate figure 1b, 1c, 1d (plots of plasmid relatedness  
% divided by chromosome relatedness). Figure 1a was made on BioRender.

close all
clearvars
clc

% default parameter values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcR=0:0.001:1; % defines the range of cooperative plasmid frequency values to be used on the x axis.
pd=0; % defector plasmid frequency value
N=4; % number of founder cells
beta=0.6; % plasmid transfer probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure 1b (vary N) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on % allows us to plot multiple lines on the same graph

% the following loop generates multiple lines for different N values.
for N=[2 4 6 8 ]

% evaluate the plasmid relatedness function for our chosen parameter
% values.
R =((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - (pcR + pcR.^2.*(N - 1) - (2.*beta.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N - (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2 + (beta.^2.*pcR.*(N - 1).*(N - 2).*(pcR + pd - 1).^2.*(N.*pcR - 3.*pcR + 1))./N.^2)./N)./((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - pcR + (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2);

% plot plasmid relatedness divided by chromosome relatedness (chromosome 
% relatedness is given by 1/N).
plot(pcR,R ./ (1/N),'LineWidth',1)

end

hold off

% formatting
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
ylim([1 2.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure 1c (vary beta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=4; % set the number of founder cells back to the default value.

figure % new figure panel

hold on % allows us to plot multiple lines on the same graph

% the following loop generates multiple lines for different beta values.
for beta=[0.2 0.4 0.6 0.8]

% evaluate the plasmid relatedness function for our chosen parameter
% values.
R =((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - (pcR + pcR.^2.*(N - 1) - (2.*beta.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N - (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2 + (beta.^2.*pcR.*(N - 1).*(N - 2).*(pcR + pd - 1).^2.*(N.*pcR - 3.*pcR + 1))./N.^2)./N)./((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - pcR + (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2);

% plot plasmid relatedness divided by chromosome relatedness (chromosome 
% relatedness is given by 1/N).
plot(pcR,R ./ (1/N),'LineWidth',1)

end

hold off

% formatting
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
ylim([1 2.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Figure 1d (vary defector plasmid freqeuncy) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta=0.6; % set theplasmid transfer probability back to the default value.

figure % new figure panel

hold on % allows us to plot multiple lines on the same graph

% the following loop generates multiple lines for different defector 
% plasmid frequencies.
for pd=[0 0.2 0.4 0.6]

% evaluate the plasmid relatedness function for our chosen parameter
% values.
R =((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - (pcR + pcR.^2.*(N - 1) - (2.*beta.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N - (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2 + (beta.^2.*pcR.*(N - 1).*(N - 2).*(pcR + pd - 1).^2.*(N.*pcR - 3.*pcR + 1))./N.^2)./N)./((pcR - (beta.*pcR.*(N - 1).*(pcR + pd - 1))./N).^2 - pcR + (beta.^2.*pcR.*(N - 1).*(pcR + pd - 1).*(N.*pcR - 2.*pcR + 1))./N.^2);

% plot plasmid relatedness divided by chromosome relatedness (chromosome 
% relatedness is given by 1/N).
plot(pcR(pcR<=(1-pd)),R(pcR<=(1-pd))./ (1/N),'LineWidth',1	)

end

hold off

% formatting
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
ylim([1 2.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

