
%Task5
%generates new shapes and returns them as [128 2] matrix.
function generatedShape = generateShape(b,sortedEigenVects,meanshape)

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

end
