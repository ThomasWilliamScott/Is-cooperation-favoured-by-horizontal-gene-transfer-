% This script tests the relatedness equation against simulated data. Our
% approach is to: (1) generate a simulated dataset for different model 
% parameter values (i.e. different values of pD, pC, beta and N), by 
% drawing from binomial distributions; (2) calculate relatedness for these 
% simulated datasets; (3) calculate the difference between the simulated 
% relatedness values and the predicted relatedness values based on our 
% algebraic relatedness expression, and show that this difference is 
% approximately zero across the whole parameter space. The data generated
% from this script is saved as
% "Compare_simulated_and_expected_relatedness.mat" in the GitHub
% repository. The results matrix "diff" shows that expected and simulated
% relatedness are equal.

clearvars
clc

% PARAMETER SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We define arrays that sweep through the full model parameter space.

pR        = 0.2:0.2:1; % Total plasmid frequency
plascoopR = 0.2:0.2:1; % Cooperator plasmid frequency
betaR     = 0.2:0.2:1; % Plasmid transfer probability
NR        = 2:2:10;    % Number of founders

% Number of simulated populations (needs to be large, otherwise our 
% relatedness expression would not be expected to equal our simulated data,
% owing to stochasticity. Note that the difference between simulated and 
% expected relatedness can be made even smaller by runing this script with
% a higher num value.)
num=100000; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATE SIMULATED DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty results matrix to be populated. This matrix will be populated with
% the difference between simulated and expected relatedness values. It is 4
% dimensional, meaning it records results for different pc, pd, N and beta
% values.
diff = NaN(length(pR),length(plascoopR),length(betaR),length(NR)) ; 

% This loops across p values (plasmid frequencies)
for curp = 1:length(pR)

    p = pR(curp);

% This loops across plascoop values (plasmid cooperativeness)    
for curplascoop = 1:length(plascoopR)

    plascoop = plascoopR(curplascoop);

    pc=p*plascoop;
    pd=p*(1-plascoop);

% This loops across beta values (plasmid transfer rates)        
for curbeta = 1:length(betaR)

    beta = betaR(curbeta);

% This loops across N values (number of founders)        
for curN = 1:length(NR)

    N = NR(curN);

AnyPlas=NaN(num,N); % Empty matrix to be populated with simulated data
PropCoop=NaN(num,N); % Empty matrix to be populated with simulated data
IPrime=NaN(num,N); % Empty matrix to be populated with simulated data

% We are looping over x (ranging from 1 to num), so that we get data for a 
% large number of simulated populations. We have set num to be large, 
% because stochatstic deviations from our relatedness expression would be 
% expected for lower num.

for x = 1:num

% First we populate two arrays, AnyPlas(x,:) & PropCoop(x,:), by drawing 
% from binomial distributions. This allows us to obtain simulated data, 
% which we can go on to use to verify our relatedness expression. 

% The 'AnyPlas(x,:)' array gives an array of N entries, one for each member
% of the deme. Each entry is either 1 (meaning the individual has a 
% plasmid) or 0 (meaning the individual lacks a plasmid). It is not
% specified in this array whether the plasmid is cooperative or defective.
AnyPlas(x,:) = (binornd(1,pc+pd,[1,N])); 


% The 'PropCoop(x,:)' array gives an array of N entries, one for each 
% member of the deme. Each entry is either 1 (meaning that, if the 
% individual happens to have a plasmid, it is cooperative) or 0 (meaning 
% that, if the individual happens to have a plasmid, it is defective).
PropCoop(x,:) = (binornd(1,pc/(pc+pd),[1,N]));

end

% Next, we generate our simulated I(x,:) array by multiplying our 
% AnyPlas(x,:) and PropCoop(x,:) arrays. The I(x,:) array gives an array of
% N entries, one for each member of the deme. Each entry is either 1 
% (meaning the individual has a cooperative plasmid) or 0 (meaning the 
% individual lacks a cooperative plasmid).
I = AnyPlas.*PropCoop;

% Next, we generate our simulated J(x,:) array by multiplying our 
% AnyPlas(x,:) and 1-PropCoop(x,:) arrays. The J(x,:) array gives an array 
% of N entries, one for each member of the deme. Each entry is either 1 
% (meaning the individual has a defector plasmid) or 0 (meaning the 
% individual lacks a defector plasmid).
J = AnyPlas.*(1-PropCoop);

% Next, we generate our simulated S score by summing all entries in the
% I(x,:) array. The S score gives the total number of cooperative plasmids
% in the deme.
S = sum(I,2);

% Next, we generate our simulated I1Prime score by combining the
% previously derived arrays in the following way. The I1Prime score is one
% if the focal individual (i.e. individual who occupies the first entry in 
% the simulated arrays) has a cooperative plasmid after horizontal gene
% transfer has taken place, and I1Prime is zero if the focal individual
% lacks a cooperative plasmid after horizontal gene transfer has taken
% place.
I1Prime = I(:,1)+(1-I(:,1)-J(:,1))*(beta/N).*(S-I(:,1));


% Next, we define a simulated array IPrime which is an array 
% of N entries, one for each member of the deme. Each entry is either 1 
% (meaning the individual has a cooperative plasmid after horizontal gene 
% transfer) or 0 (meaning the individual lacks a cooperative plasmid after 
% horizontal gene transfer).
for y=1:N
IPrime(:,y) = I(:,y)+(1-I(:,y)-J(:,y))*(beta/N).*(S-I(:,y));
end

% The simulated relatedness value is therefore calculated using the
% following formula, which is a general formula for relatedness, derived in
% Grafen (1985) and elsewhere. We need to extract the actual covariance
% (which we call simulated_rel) from the covariance matrix (which we call 
% simulated_rel_mat). 
simulated_rel_mat = cov(mean(IPrime,2),I1Prime) ./ var(I1Prime);
simulated_rel = simulated_rel_mat(1,2);

% Our algebraic relatedness expression is given below.
algebraic_rel = ((pc - (beta*pc*(N - 1)*(pc + pd - 1))/N)^2 - (pc + pc^2*(N - 1) - (2*beta*pc*(N - 1)*(pc + pd - 1)*(N*pc - 2*pc + 1))/N - (beta^2*pc*(N - 1)*(pc + pd - 1)*(N*pc - 2*pc + 1))/N^2 + (beta^2*pc*(N - 1)*(N - 2)*(pc + pd - 1)^2*(N*pc - 3*pc + 1))/N^2)/N)/((pc - (beta*pc*(N - 1)*(pc + pd - 1))/N)^2 - pc + (beta^2*pc*(N - 1)*(pc + pd - 1)*(N*pc - 2*pc + 1))/N^2);

% We record the difference between the simulated and algebraic relatedness
% below. A difference of approximately zero implies that our alegbraic 
% relatedness expression is correct.

diff(curp,curplascoop,curbeta,curN) = algebraic_rel - simulated_rel;

end
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We saved the output of this script using the following line of code 
% (commented out):
        %save("Comparison_of_simulated_and_expected_relatedness.mat")

% By looking at the diff matrix, it can be seen that the highest recorded
% difference between expected and simulated relatedness (given by 
% max(max(max(max((abs(diff))))))) is 0.0052. This indicates a good 
% correspondence between simulated and expected relatedness. Note that this
% difference can be made even lower by increasing the 'num' parameter.
