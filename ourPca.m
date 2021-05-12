
%computes pca
function [geordneteEigenVects, geordneteEigenVals] = ourPca(data)
    
    data = data - mean(data,2);

    [geordneteEigenVals,i] = sort(eig(ourCov(data)),'descend');

    [vect,val] = eig(ourCov(data));
    geordneteEigenVects = vect(:,i);
    return
end