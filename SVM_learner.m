function model = SVM_learner(hs,traindata, trainlabels, C, kerneltype, r)
    % INPUT : 
    % traindata   - m X n matrix, where m is the number of training points
    % trainlabels - m X 1 vector of training labels for the training data
    % C           - SVM regularization parameter (positive real number)
    % kerneltype  - one of strings 'linear', 'poly', 'rbf'
    %               corresponding to linear, polynomial, and RBF kernels
    %               respectively.
    % r           - integer parameter indicating the degree of the
    %               polynomial for polynomial kernel, or the width
    %               parameter for the RBF kernel; not used in the case of
    %               linear kerne and can be set to a default value.
    
    % OUTPUT
    % returns the structure 'model' which has the following fields, in
    % addition to the training data/parameters.(You can choose to add more
    % fields to this structure needed for your implementation)
    
    
    % 	alphas      	- m X 1 vector of support vector coefficients
    % 	b           	- SVM bias term
    % 	objective   	- optimal objective value of the SVM solver
    % 	support_vectors - the subset of training data, which are the support vectors
    %   w               - the classifier vector
    
    % Default code below. Fill in the code for solving the
    % SVM dual optimization problem using quadprog function
%     disp(kerneltype);
%      disp(r);
    b = 0;
    objective = 0;
    s = size(traindata, 1);
    alphas = repmat(0, s, 1);
    model.kerneltype = kerneltype;
    model.r = r;
    model.C = C;
    model.traindata = traindata;
    mode.trainlabels = trainlabels;
    support_vectors = [];
    %k = repmat(0, s, s);
    %for i = 1 : s
     %   for j = 1 : s
            k = compute_kernel(traindata, traindata, kerneltype, r);
        %end
    %end
    H = k.*(trainlabels*trainlabels');
    Aeq = trainlabels.';
    %beq = 0;
    beq = 0;
    %disp(trainlabels);
    f = ones (s,1);
    A = [];
    b = [];
    %lb = 0;
    lb = zeros(s,1);
    ub = C*ones (s,1);
    options = optimset('Algorithm','interior-point-convex','display','off','MaxIter',10000);
    if(hs=='h')
    [alphas, model.objective] = quadprog(H,-f,A,b,Aeq,beq,lb,[],[],options);
    else
    [alphas, model.objective] = quadprog(H,-f,A,b,Aeq,beq,lb,ub,[],options);
    end
    model.objective = -1*model.objective;
    %disp(model.objective);
    w = zeros(1,size(traindata,2));
    for i= 1 : s
       w = w + alphas(i,:)*trainlabels(i,:)*traindata(i,:); 
    end
    
    model.w = w;
    %disp(support_vectors);
end
