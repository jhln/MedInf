
%calls generateShape.m and plots the result.
function plotShape(meanshape,eigenVects,b,color)
%plot mean shape
meanB = zeros(length(b),1);
generatedShape = generateShape(meanB,eigenVects,meanshape);
plot(generatedShape(1,:),generatedShape(2,:),'r');

%plot generated shape 
generatedShape = real(generateShape(b,eigenVects,meanshape));
pl = plot(generatedShape(1,:),generatedShape(2,:),color);
%pl.Color(4) = 0.2;
end
