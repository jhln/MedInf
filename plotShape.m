
function plotShape(shapes)
%plot mean shape
ones = [1; 1; 1];
[generatedShape,meanshape] = generateShape(ones,shapes);
plot(meanshape(1,:),meanshape(2,:),'r');
axis equal
hold on
b = [100; 0; 0];
[generatedShape,meanshape] = generateShape(b,shapes);
plot(generatedShape(1,:),generatedShape(2,:),'b');

b = [0; 100; 0];
[generatedShape,meanshape] = generateShape(b,shapes);
plot(generatedShape(1,:),generatedShape(2,:),'g');

b = [0; 0; 100];
[generatedShape,meanshape] = generateShape(b,shapes);
plot(generatedShape(1,:),generatedShape(2,:),'c');

hold off
end

function plotRandomShape(shapes)
%plot mean shape
ones = ones(3,1);
[generatedShape,meanshape] = generateShape(ones,shapes);
plot(meanshape(1,:),meanshape(2,:),'r');
axis equal
hold on
%plot random new shapes
for i = 1:2  
    ones = (randn(3,1));
    ones = ones/norm(ones);
    ones = ones .* 100;
    ones = [0; 100; 0];
    [generatedShape,meanshape] = generateShape(ones,shapes);
    plot(generatedShape(1,:),generatedShape(2,:),'b');
end
hold off
end