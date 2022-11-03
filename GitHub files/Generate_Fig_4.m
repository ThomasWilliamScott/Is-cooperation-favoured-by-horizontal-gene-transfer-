% This script generates the two 'single trial' panels in Figure 4. The
% remaining 'overview' panel in Figure 4 was made in powerpoint using the
% data saved in "Eq_Genotype_Freqs_N=20.mat" and plotted in heatmap form
% in Figure 2.

clearvars
close all
clc

% Specify parameter inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=500; % Generations.
B=1.435; % Benefit of cooperation.
CG=0.1; % Cost of cooperation.
CC=0.2; % Cost of plasmid carriage. 
N = 20; % Number of founders.
beta=0.95; % Plasmid transfer rate.

% specify initial genotype frequencies
x12(1) = 0.0001; % Initial frequency of genotype 12.
x11(1) = 0.0001; % Initial frequency of genotype 11.
x22(1) = 0.0001; % Initial frequency of genotype 22.
x21(1) = 0.0001; % Initial frequency of genotype 21.
x20(1) = 0.0001; % Initial frequency of genotype 20.
x10(1) = 0.9995;  % Initial frequency of genotype 10 (chromosomal defectors initially common).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for s=[0.1 0.3]; % Plasmid loss rate (we plot graphs for two different values)


% The following code iterates the population genetic recursions described
% in the Appendix. It is copied from the 
% "FINAL_Script_to_generate_equilibrium_genotype_frequencies" script.

for t=1:T % Iterate recursions over T generations.   
    
    % Using for loops, we calculate patch cooperator frequencies,
    % for members of each genotype. In other words, here, we are
    % calculating the m_{A} terms described in the Appendix.
    
    clear m
    
    for i=[0:2]  % no. of 20,21 individuals
        for j=[0:2] % no. of 22,12 individuals
            for k=[0:2] % no. of 10 individuals
                for l=[0:2] % no. of 11 individuals
                    if i+j+k+l==2 || i+j+k+l==1
                          K=i+j+k+l;
                          m(i+1,j+1,k+1,l+1)= ( N*(i+j)+ beta*j*k + (K-N)*(beta*x10(t)*(-j+x12(t)*(K-N+1)+x22(t)*(K-N+1))+x12(t)*(-(beta*k+N))-x22(t)*(beta*k+N)-N*(x20(t)+x21(t))) )/(N^2) ;
                    end
                end
            end
        end
    end
    m10= m(0+1,0+1,1+1,0+1);
    m11= m(0+1,0+1,0+1,1+1);
    m12= m(0+1,1+1,0+1,0+1);
    m20= m(1+1,0+1,0+1,0+1);
    m21= m(1+1,0+1,0+1,0+1);
    m22= m(0+1,1+1,0+1,0+1);
    
    m2210= m(0+1,1+1,1+1,0+1);
    m2211= m(0+1,1+1,0+1,1+1);
    m2212= m(0+1,2+1,0+1,0+1);
    m2220= m(1+1,1+1,0+1,0+1);
    m2221= m(1+1,1+1,0+1,0+1);
    m2222= m(0+1,2+1,0+1,0+1);
    
    m2110= m(1+1,0+1,1+1,0+1);
    m2111= m(1+1,0+1,0+1,1+1);
    m2112= m(1+1,1+1,0+1,0+1);
    m2121= m(2+1,0+1,0+1,0+1);
    m2120= m(2+1,0+1,0+1,0+1);
    
    m2010= m(1+1,0+1,1+1,0+1);
    m2011= m(1+1,0+1,0+1,1+1);
    m2012= m(1+1,1+1,0+1,0+1);
    m2020= m(2+1,0+1,0+1,0+1);

    m1210= m(0+1,1+1,1+1,0+1);
    m1211= m(0+1,1+1,0+1,1+1);
    m1212= m(0+1,2+1,0+1,0+1);

    m1110= m(0+1,0+1,1+1,1+1);
    m1111= m(0+1,0+1,0+1,2+1);
   
    m1010= m(0+1,0+1,2+1,0+1);  
    
    % Having calculated patch cooperator frequencies, we can now explicitly
    % calculate mean population fitness (equation also given in the
    % Appendix).
    
    W(t) = (x22(t)*((1/N)*(1+m22*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2222*B-CG-CC)+x21(t)*(1+m2221*B-CG-CC)+x20(t)*(1+beta)*(1+m2220*B-CG-CC)+x12(t)*(1+m2212*B-CG-CC)+x11(t)*(1+m2211*B-CG-CC)+x10(t)*(1+m2210*B-CG-CC)))+beta*((N-1)/N)*x20(t)*x12(t)*(1+m2012*B-CG-CC))+(x21(t)*((1/N)*(1+m21*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2221*B-CG-CC)+x21(t)*(1+m2121*B-CG-CC)+x20(t)*(1+beta)*(1+m2120*B-CG-CC)+x12(t)*(1+m2112*B-CG-CC)+x11(t)*(1+m2111*B-CG-CC)+x10(t)*(1+m2110*B-CG-CC)))+beta*((N-1)/N)*x20(t)*x11(t)*(1+m2011*B-CG-CC))+(x20(t)*((1/N)*(1+m20*B-CG)+((N-1)/N)*(x22(t)*(1-beta)*(1+m2220*B-CG)+x21(t)*(1-beta)*(1+m2120*B-CG)+x20(t)*(1+m2020*B-CG)+x12(t)*(1-beta)*(1+m2012*B-CG)+x11(t)*(1-beta)*(1+m2011*B-CG)+x10(t)*(1+m2010*B-CG))))+(x12(t)*((1/N)*(1+m12*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2212*B-CG-CC)+x21(t)*(1+m2112*B-CG-CC)+x20(t)*(1+m2012*B-CG-CC)+x12(t)*(1+m1212*B-CG-CC)+x11(t)*(1+m1211*B-CG-CC)+x10(t)*(1+beta)*(1+m1210*B-CG-CC)))+beta*((N-1)/N)*x22(t)*x10(t)*(1+m2210*B-CG-CC))+(x11(t)*((1/N)*(1+m11*B-CC)+((N-1)/N)*(x22(t)*(1+m2211*B-CC)+x21(t)*(1+m2111*B-CC)+x20(t)*(1+m2011*B-CC)+x12(t)*(1+m1211*B-CC)+x11(t)*(1+m1111*B-CC)+x10(t)*(1+beta)*(1+m1110*B-CC)))+beta*((N-1)/N)*x21(t)*x10(t)*(1+m2110*B-CC))+(x10(t)*((1/N)*(1+m10*B)   +((N-1)/N)*(x22(t)*(1-beta)*(1+m2210*B)   +x21(t)*(1-beta)*(1+m2110*B)   +x20(t)*(1+m2010*B)   +x12(t)*(1-beta)*(1+m1210*B)   +x11(t)*(1-beta)*(1+m1110*B)   +x10(t)*(1+m1010*B))));

    % Having calculated patch cooperator frequencies and mean population
    % fitness, we can now explicitly calculate the frequency of each
    % genotype after selection (and before plasmid loss).
    
    x22M = (x22(t)*((1/N)*(1+m22*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2222*B-CG-CC)+x21(t)*(1+m2221*B-CG-CC)+x20(t)*(1+beta)*(1+m2220*B-CG-CC)+x12(t)*(1+m2212*B-CG-CC)+x11(t)*(1+m2211*B-CG-CC)+x10(t)*(1+m2210*B-CG-CC)))+beta*((N-1)/N)*x20(t)*x12(t)*(1+m2012*B-CG-CC))/W(t); 
    x21M = (x21(t)*((1/N)*(1+m21*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2221*B-CG-CC)+x21(t)*(1+m2121*B-CG-CC)+x20(t)*(1+beta)*(1+m2120*B-CG-CC)+x12(t)*(1+m2112*B-CG-CC)+x11(t)*(1+m2111*B-CG-CC)+x10(t)*(1+m2110*B-CG-CC)))+beta*((N-1)/N)*x20(t)*x11(t)*(1+m2011*B-CG-CC))/W(t);  
    x20M = (x20(t)*((1/N)*(1+m20*B-CG)+((N-1)/N)*(x22(t)*(1-beta)*(1+m2220*B-CG)+x21(t)*(1-beta)*(1+m2120*B-CG)+x20(t)*(1+m2020*B-CG)+x12(t)*(1-beta)*(1+m2012*B-CG)+x11(t)*(1-beta)*(1+m2011*B-CG)+x10(t)*(1+m2010*B-CG))))/W(t);
    x12M = (x12(t)*((1/N)*(1+m12*B-CG-CC)+((N-1)/N)*(x22(t)*(1+m2212*B-CG-CC)+x21(t)*(1+m2112*B-CG-CC)+x20(t)*(1+m2012*B-CG-CC)+x12(t)*(1+m1212*B-CG-CC)+x11(t)*(1+m1211*B-CG-CC)+x10(t)*(1+beta)*(1+m1210*B-CG-CC)))+beta*((N-1)/N)*x22(t)*x10(t)*(1+m2210*B-CG-CC))/W(t); 
    x11M = (x11(t)*((1/N)*(1+m11*B-CC)+((N-1)/N)*(x22(t)*(1+m2211*B-CC)+x21(t)*(1+m2111*B-CC)+x20(t)*(1+m2011*B-CC)+x12(t)*(1+m1211*B-CC)+x11(t)*(1+m1111*B-CC)+x10(t)*(1+beta)*(1+m1110*B-CC)))+beta*((N-1)/N)*x21(t)*x10(t)*(1+m2110*B-CC))/W(t);
    x10M = (x10(t)*((1/N)*(1+m10*B)   +((N-1)/N)*(x22(t)*(1-beta)*(1+m2210*B)   +x21(t)*(1-beta)*(1+m2110*B)   +x20(t)*(1+m2010*B)   +x12(t)*(1-beta)*(1+m1210*B)   +x11(t)*(1-beta)*(1+m1110*B)   +x10(t)*(1+m1010*B))))/W(t);
      
    % And finally, we can calculate the freqeuncy of each genotype after
    % plasmid loss. This gives the frequency of each genotype at the start 
    % of the next generation, which is why we enter these values into the 
    % (t+1) position in the genotype frequency arrays.
    x22(t+1) = x22M*(1-s);
    x21(t+1) = x21M*(1-s);
    x20(t+1) = x20M + x21M*s + x22M*s;
    x12(t+1) = x12M*(1-s);
    x11(t+1) = x11M*(1-s);
    x10(t+1) = x10M + x12M*s + x11M*s;
    
end

% This records the plasmid frequency over time.
res_plas = x22+x21+x12+x11 ; 

% This records the cooperator frequency over time.
res_coop = x21+x22+x12+x20 ; 

% This records plasmid relatedness over time.
res_rel = (((x12 + x22) - (beta.*(x12 + x22).*(N - 1).*((x12 + x22) + (x11 + x21) - 1))./N).^2 - ((x12 + x22) + (x12 + x22).^2.*(N - 1) - (2.*beta.*(x12 + x22).*(N - 1).*((x12 + x22) + (x11 + x21) - 1).*(N.*(x12 + x22) - 2.*(x12 + x22) + 1))./N - (beta.^2.*(x12 + x22).*(N - 1).*((x12 + x22) + (x11 + x21) - 1).*(N.*(x12 + x22) - 2.*(x12 + x22) + 1))./N.^2 + (beta.^2.*(x12 + x22).*(N - 1).*(N - 2).*((x12 + x22) + (x11 + x21) - 1).^2.*(N.*(x12 + x22) - 3.*(x12 + x22) + 1))./N.^2)./N)./(((x12 + x22) - (beta.*(x12 + x22).*(N - 1).*((x12 + x22) + (x11 + x21) - 1))./N).^2 - (x12 + x22) + (beta.^2.*(x12 + x22).*(N - 1).*((x12 + x22) + (x11 + x21) - 1).*(N.*(x12 + x22) - 2.*(x12 + x22) + 1))./N.^2);

% This records plasmid cooperativeness over time.
res_plas_coop = (x22+x12) ./ res_plas;

% This sets plasmid cooperativeness to undefined when the plasmid is absent.
res_plas_coop(res_plas<0.00001)=NaN;

% This sets plasmid relatedness to undefined when the plasmid is absent.
res_rel(res_plas<0.00001)=NaN;


figure % new figure panel

hold on % plot multiple lines on same figure

% plots freuqency of the 12 genotype
plot(1:T+1,x12,'LineWidth',2)

% plots freuqency of the 11 genotype
plot(1:T+1,x11,'LineWidth',2)

% plots plasmid relatedness
plot(1:T+1,res_rel,'LineWidth',2,'LineStyle','--','Color',[0.3 0.3 0.3])
hold off

% formatting 
xlim([1 T+1])
ylim([0 1])
box off
set(gcf,'color','white')
set(gca,'fontsize',16)
yticks([0:0.2:1])

end