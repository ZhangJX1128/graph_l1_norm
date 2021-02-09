#graph_l1_norm

Prediction of lncRNA-disease associations based on robust multi-label learning

Method Description:
we present a novel robust computational framework for lncRNA-disease association prediction by combining ℓ1-norm graph with multi-label learning. Specifically, we first construct a set of similarity matrices for lncRNAs and diseases using known associations. Then, both lncRNA and disease similarity matrices are adaptively re-weighted to enhance the robustness via the ℓ1-norm graph. Lastly, the association matrix is updated with a graph-based multi-label learning framework to uncover the underlying consistency between the lncRNA space and the disease space.

Requirements:
MCGLLDA was developed in MATLAB 2016b environment, but it should be working in all MATLAB versions.

Usage:
1.Dataset
lncSim.mat and disSim_Jaccard.mat store lncRNA similarity matrix and disease similarity matrix, respectively; interMatrix.mat stores known lncRNA-disease association information;lncRNA_Name.txt and diseases_Name.txt store lncRNA ids and disease ids, respectively.
2.Code
To run our method, simply open the "caseStudy.m" script in matlab programming environment and press "Run" button,you will get the predict potential lncRNA-disease associations.

Contact:
For any questions regarding our work, please feel free to contact us: alcs417@sdnu.edu.cn.

