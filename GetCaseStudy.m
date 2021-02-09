function [F ] = GetCaseStudy( lncR_sim_matrix, disease_sim_matrix, lncR_disease_matrix )
        SL = lncR_sim_matrix;
        SD = disease_sim_matrix;
        Y_temp = lncR_disease_matrix;
        
        [KD1, KL1] = GaussianKernel(Y_temp', 1, 1);
        [KD2, KL2] = Cosine(Y_temp');
        
        SL = (SL + KL1 + KL2) / 3; 
        SD = (SD + KD1 + KD2) / 3; 
        
        F = graph_l1_norm(SL, SD, Y_temp, 1);
%         save F;

end

