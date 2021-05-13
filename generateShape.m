
function [generatedShape,meanShape2] = generateShape(b,data)

formattedData = format(data);
meanShape = mean(formattedData,1);
formattedData = formattedData - meanShape;
%mymax = max(formattedData2,[],'all');
%formattedData2 = formattedData ./ max(abs(formattedData),[],'all');
%meanShape = meanShape ./ (max(abs(meanShape),[],'all'));


[sortedEigVectsShapes, sortedEigValsShapes] = ourPca(transpose(formattedData));
numberEigenVectors = length(b);

temp = (sortedEigVectsShapes(:,1:numberEigenVectors));
temp2 = temp * b;
gs = temp2 + transpose(meanShape);
generatedShape = zeros(2,128);
meanShape2 = zeros(2,128);
index = 1;
for i = 1:128
generatedShape(1,i) = gs(index,1);
meanShape2(1,i) = meanShape(1,index);
index = index + 1;
generatedShape(2,i) = gs(index,1);
meanShape2(2,i) = meanShape(1,index);
index = index + 1;
end



end

%formate bone shapes into 14x256 matrix (one row = x1 y1 x2 y2 ....)
function formattedShapes = format(data)
formattedShapes = zeros(14,256);
for j = 1:14
temp = zeros(1,256);
index = 1;
for i = 1:128
    temp(1,index) = data(i,1,j);
    index = index + 1;
    temp(1,index) = data(i,2,j);
    index = index + 1;
end
formattedShapes(j,:) = temp;
end
end
