function geordneteEigenVects = pca(data)
    %[v,d] = eig(data)
    
    [eigenValues,i] = sort(eig(ourCov2(data)),'descend')
    [vect,val] = eig(ourCov2(data))

    geordneteEigenVects = vect(:,i)
end