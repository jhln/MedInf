
%HXMzNnvuht7WBgJ

%data3D = load('daten3d.mat');

%min_and_max = @(x) [min(x), max(x)];
%m = min_and_max([3 4 1 6 2])


%predefined covarianceMatrix Functions
%m2D = cov(data2D.data1(1,:),data2D.data1(2,:))
%m3D = cov(data3D.data)


%mat11 = covFunc(data.data1(:,1),data.data1(:,1));
%mat12 = covFunc(data.data1(:,1),data.data1(:,2));
%mat21 = covFunc(data.data1(:,2),data.data1(:,1));
%mat22 = covFunc(data.data1(:,2),data.data1(:,2));

%ownMat2D = covMat(data2D.data1)
%ownMat3D = covMat(data3D.data)

%[r,c] = size(dataMatrix)


% self-programmed covarianceMatrix with variable dimensions
function covMatrix = ourCov2(dataMatrix)

    [r,c] = size(dataMatrix);
    
    covMatrix = ones(r);
    
    for ix = 1:r
        for iy = 1:r
            covMatrix(ix,iy) = covFunc(dataMatrix(ix,:),dataMatrix(iy,:));
        end
    end

    %mat11 = covFunc(dataMatrix(:,1),dataMatrix(:,1));
    %mat12 = covFunc(dataMatrix(:,1),dataMatrix(:,2));
    %mat21 = covFunc(dataMatrix(:,2),dataMatrix(:,1));
    %mat22 = covFunc(dataMatrix(:,2),dataMatrix(:,2));

    %covMatrix = [ mat11 mat12; mat21 mat22]
    return
end
    

% helperFunction
function result = covFunc(column1, column2)
    %data = load('/home/johannes/uni/21ss/MedInf/daten.mat');
    %m = cov(data.data1(:,1),data.data1(:,2));
    %disp(m);
    
    weight = 1/(length(column1) - 1 );
    
    mean1 = mean(column1);
    mean2 = mean(column2);
    
    result = weight * (sum((column1 - mean1) .* (column2 - mean2)));
    return
end


