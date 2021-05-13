
%computes pca
function [geordneteEigenVects, geordneteEigenVals] = ourPca(data)
    
    [geordneteEigenVals,i] = sort(eig(ourCov(data)),'descend');

    [vect,val] = eig(ourCov(data));
    geordneteEigenVects = vect(:,i);
    return
end