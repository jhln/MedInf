
% self-programmed covarianceMatrix with variable dimensions
function covMatrix = ourCov(dataMatrix)

    [r,c] = size(dataMatrix);
    
    covMatrix = ones(r);
    
    for ix = 1:r
        for iy = 1:r
            covMatrix(ix,iy) = covFunc(dataMatrix(ix,:),dataMatrix(iy,:));
        end
    end
    return
end
    

% helperFunction
function result = covFunc(column1, column2)
    
    weight = 1/(length(column1) - 1 );
    
    mean1 = mean(column1);
    mean2 = mean(column2);
    %mean1 = 0;
    %mean2 = 0;
    
    result = weight * (sum((column1-mean1) .* (column2-mean2)));
    return
end


