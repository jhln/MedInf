
%Task5
%generates new shapes and returns them as [128 2] matrix.
function generatedShape = generateShape(b,s,r,x,y,sortedEigenVects,meanshape)

    numberEigenVectors = length(b);

    temp = (sortedEigenVects(:,1:numberEigenVectors));
    temp2 = temp * b;
    gs = temp2 + transpose(meanshape);
    generatedShape = zeros(2,128);
    index = 1;
    for i = 1:128
        generatedShape(1,i) = gs(index,1);
        index = index + 1;
        generatedShape(2,i) = gs(index,1);
        index = index + 1;
    end
    sz = size(generatedShape,1);
    sm = scalingMatrix(s,[sz sz]);
    generatedShape = sm * generatedShape;
    rm = rotationMatrix(deg2rad(r));
    generatedShape = rm * generatedShape;
    generatedShape(1,:) = generatedShape(1,:) + x;
    generatedShape(2,:) = generatedShape(2,:) + y;
end

function sm = scalingMatrix(s,size)
    sm = eye(size);
    sm = sm .* s;
    %sm(size, size) = 1;
end

function rm = rotationMatrix(r)
    rm = ones(2,2);
    rm(1,1) = cos(r);
    rm(1,2) = -sin(r);
    rm(2,1) = sin(r);
    rm(2,2) = cos(r);
end
