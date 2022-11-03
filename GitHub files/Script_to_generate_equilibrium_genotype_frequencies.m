% This script solves the theoretical population genetic model described in
% the Appendix. More specifically, in the Appendix, we derived a series
% of recursion equations, which describe how genotype frequencies change
% over a single generation. This means that we can iterate these recursion
% equations to track the evolutionary process (see how genotype frequencies
% change over time) and find the equilibrium state of the population. This
% script was written for the purpose of iterating our recursion equations.
% Iterating this program generates six outputs: res_x22, res_x21, res_x20, 
% res_x12, res_x11, res_x10. These outputs show the equilibrium frequency 
% of each genotype (respectively, 22, 21, 20, 12, 11, 10). All outputs are 
% matrices, giving equilibrium frequencies across a range of beta (plasmid 
% transfer rate) and s (plasmid loss) values. These outputs form the 
% basis of all 'heatmap' figures (Fig. 2 / 4, Supp. Fig. 1).

clearvars
close all
clc

% Specify parameter inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=100000; % Generations.
B=1.435; % Benefit of cooperation.
CG=0.1; % Cost of cooperation.
CC=0.2; % Cost of plasmid carriage. 
betaR = [0.0:0.018:1]; % Plasmid transfer rate (we generate data for a range of rates).
N = 20; % Number of founders (this is the parameter value used for Fig. 2 / 4; it was changed to N=15 for Supp. Fig. 1).

% Set initial genotype frequencies so that genotype 10 (defector, no 
% plasmid) is initially at 0.9 frequency, and the remaining 0.1 frequency 
% is distributed amongst the remaining genotypes.
x12(1) = rand; 
x11(1) = rand; 
x22(1) = rand; 
x21(1) = rand; 
x20(1) = rand; 
x10(1) = 0.9;  % Initial frequency of genotype 10 (chromosomal defectors initially common).
tot = 10*(sum(x11(1)+x12(1)+x20(1)+x21(1)+x22(1)));
x12(1) = x12(1) / tot; % Initial frequency of genotype 12.
x11(1) = x11(1) / tot; % Initial frequency of genotype 11.
x22(1) = x22(1) / tot; % Initial frequency of genotype 22.
x21(1) = x21(1) / tot; % Initial frequency of genotype 21.
x20(1) = x20(1) / tot; % Initial frequency of genotype 20.

sR=0:0.01:0.5; % range of plasmid loss values to collect data for.

% define 'NaN' matrices to populate with results.
res_x22 = NaN(length(betaR),length(sR));
res_x21 = NaN(length(betaR),length(sR));
res_x20 = NaN(length(betaR),length(sR));
res_x12 = NaN(length(betaR),length(sR));
res_x11 = NaN(length(betaR),length(sR));
res_x10 = NaN(length(betaR),length(sR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Run the population genetic model for different s and beta values %%%%%%%%

for cur_s = 1:numel(sR) % Loop over range of s values.
    
     s = sR(cur_s);
    
 for cur_beta = 1:numel(betaR) % Loop over range of beta values.
     
     beta= betaR(cur_beta) ;
     
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
    % Appendix)
    
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

% The following matrices record the equilibrium genotype frequencies, for 
% different beta and s values.
    res_x22(cur_beta,cur_s) = x22(T);
    res_x21(cur_beta,cur_s) = x21(T);
    res_x20(cur_beta,cur_s) = x20(T);
    res_x12(cur_beta,cur_s) = x12(T);
    res_x11(cur_beta,cur_s) = x11(T);
    res_x10(cur_beta,cur_s) = x10(T);
  
end
end

% The following line of code (which is commented out) was used to save the
% results. The saved datasets can be found in the GitHub repository. They
% are loaded in other scripts to generate figures.

       % save("Eq_Genotype_Freqs_N="+N+".mat")

