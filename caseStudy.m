%  Prediction of lncRNA-disease associations based on inductive matrix completion
%  

%% configuration
addpath('dataset');
 
%% load data
%%
%   interMatrix.mat: an n*m association matrix between lncRNAs and diseases, n is 
%    the number of lncRNAs, and m is the number of diseases
%   lncSim.mat: an n*n sequence similarity matrix of lncRNAs
%   disSim_Jaccard.mat: an m*m similarity matrix of disease

   load(sprintf('./datasets/%s/disSim_Jaccard.mat', dataPath));
   load(sprintf('./datasets/%s/interMatrix.mat', dataPath));
   load(sprintf('./datasets/%s/lncSim.mat', dataPath));


% lncR_sim_matrix = lncSim;
% disease_sim_matrix = disSim_Jaccard;
% lncR_disease_matrix = interMatrix; 



prediction = GetCaseStudy(lncR_sim_matrix, disease_sim_matrix, lncR_disease_matrix);

  save('output/caseStudy.mat','prediction');

 
