
%calls generateShape.m and plots the result.
function plotShape(meanshape,eigenVects,b,s,r,x,y,color)
    %plot mean shape
    meanB = zeros(length(b),1);
    generatedShape = generateShape(meanB,s,r,x,y,eigenVects,meanshape);
    plot(generatedShape(1,:),generatedShape(2,:),'r');

    %plot generated shape 
    generatedShape = real(generateShape(b,s,r,x,y,eigenVects,meanshape));
    pl = plot(generatedShape(1,:),generatedShape(2,:),color);
    %pl.Color(4) = 0.2;
end
