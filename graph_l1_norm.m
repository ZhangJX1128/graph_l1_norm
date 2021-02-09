
function [F, obj] = graph_l1_norm(WL, WD, Y, ~)
    % Semi-supervised Learning using L21-norm minimization
    % min_F  0.5*\sum_ij A_ij*||Fi-Fj||_2 + tr(F-Y)'*U*(F-Y)
    
    % WL: lncRNA similarity matrix, n*n
    % WD: disease similarity matrix, m*m
    % Y: known lncRNA disease association matrix, n*m 
    % F: predicted association matrix, n*m matrix
    % obj: objective values during the iterations
    
    ite_num = 30;

    obj = zeros(ite_num,1);  
    n = size(WL, 1);
    %u = ones(n, 1);
     u = zeros(n, 1);
   
    if nargin == 4
        y = sum(Y, 2);
        u(y >= 1) =100;
    end

    U = diag(u);

    % unnormalized Laplacian matrix
    DL = diag(sum(WL, 2));
    DD = diag(sum(WD, 2));
    LL = DL - (WL + WL') / 2;
    LD = DD - (WD + WD') / 2;
      
    % use the normalized Laplacian matrix
%     WL0 = WL - diag(diag(WL));
%     WL1 = (WL0 + WL0') / 2;
%     DL = diag(sum(WL1));
%     idL = eye(size(WL1, 1));
%     WL2 = diag(1 ./ diag(sqrtm(DL)));
%     LL = idL - WL2 * WL1 * WL2;
%     
%     WD0 = WD - diag(diag(WD));
%     WD1 = (WD0 + WD0') / 2;
%     DD = diag(sum(WD1));
%     idD = eye(size(WD1, 1));
%     WD2 = diag(1 ./ diag(sqrtm(DD)));
%     LD = idD - WD2 * WD1 * WD2;
        
    F = sylvester(LL + U, LD, 0.5 * U * Y);
%     F = sylvester(2*LL + U, 2*LD, U * Y);
%     F = sylvester(LL + U, LD, U * Y);
    distQL = sqrt(abs(L2_distance_subfun(F', F')) + eps);
    distQD = sqrt(abs(L2_distance_subfun(F, F)) + eps);

    for iter = 1:ite_num
        WLL1 = WL ./ (2 * distQL);
        WDD1 = WD ./ (2 * distQD);
        
        % use the unnormalized version
        WLL1 = (WLL1 + WLL1') / 2;
        WDD1 = (WDD1 + WDD1') / 2;
        DL1 = diag(sum(WLL1, 2));
        DD1 = diag(sum(WDD1, 2));
        LL1 = DL1 - WLL1;
        LD1 = DD1 - WDD1;
        
        % use the normalized Laplacian matrix
%         WL1 = (WLL1 + WLL1') / 2;
%         DL = diag(sum(WL1));
%         idL = eye(size(WL1, 1));
%         WL2 = diag(1 ./ diag(sqrtm(DL)));
%         LL1 = idL - WL2 * WL1 * WL2;
% 
%         WD1 = (WDD1 + WDD1') / 2;
%         DD = diag(sum(WD1));
%         idD = eye(size(WD1, 1));
%         WD2 = diag(1 ./ diag(sqrtm(DD)));
%         LD1 = idD - WD2 * WD1 * WD2;
        
        F = sylvester(LL1 + U, LD1, 0.5 * U * Y);                                                                                                                                                
%           F = sylvester(2*LL1 + U, 2*LD1, U * Y);
%           F = sylvester(LL1 + U, LD1, U * Y);
        distQL = sqrt(abs(L2_distance_subfun(F', F')) + eps);
        distQD = sqrt(abs(L2_distance_subfun(F, F)) + eps);
        
        obj(iter) = 0.5 * sum(sum((distQL .* WL)))+ 0.5 * sum(sum((distQD .* WD)))+ trace((F-Y)' * U * (F-Y))
%         obj(iter) = sum(sum((distQL .* WL))) 
%                         + sum(sum((distQD .* WD))) 
%                         + trace((F-Y)' * U * (F-Y)); 

    end

end



