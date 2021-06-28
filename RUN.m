
clear all

%data load
data2D = load('daten.mat');
data3D = load('daten3d.mat');
dataS = load('shapes.mat');

%data to use for the run
activeData = data2D.data3;


%
%   Aufgabe 1
%
%activeData = activeData - mean(activeData,2);
covarianzMatrix = ourCov(activeData);

%iterate over all data matrices in e.g. data2D
%{
fn = fieldnames(data2D);
c = numel(fn);

tiledlayout(c,1)
for i = 1:c
    d = data2D.(fn{i});
    d = d - mean(d,2);
    
    subtract the mean of the data 
    d = activeData - mean(activeData,2);

    nexttile
    plot(d(1,:), d(2,:),'x')
    axis equal

    title('data')
    cov(d(1,:),d(2,:))
    ourCov(d)
end
%}

figure
plot(activeData(1,:), activeData(2,:),'x')
axis equal
title('Task 1')


%
%   Aufgabe 2
%

[eigenVects, eigenVals] = ourPca(activeData);
hauptVect = eigenVects(:,1);
plot2DPCA(transpose(activeData),mean(activeData,2),0,eigenVects,diag(eigenVals),1,0)
axis equal
title('Task 2')


%
%   Aufgabe 3
%
%shift data around mean done in task1. 

[eigenVects, eigenVals] = ourPca(activeData);
hauptVect = eigenVects(:,1);
%helperProjection = projectOnHauptVect(activeData, hauptVect);

%flip the eigenvector sign to get representative data 
projectedData = transpose(activeData) * -hauptVect;

reconstructedData = -hauptVect * transpose(projectedData(:,1)) + mean(activeData,2);
%additionally ploting the projected data (+)
plot2DPCA(transpose(activeData),mean(activeData,2),(reconstructedData)',eigenVects,diag(eigenVals),1,1);
%plot(helperProjection(1,:),helperProjection(2,:),'o')
plot(projectedData(:,1),zeros(length(projectedData(:,1)),2),'+')
title('Task 3a')

meanSquaredError = mean((activeData - reconstructedData).^2,'all');

%projection onto sidevector
%{
nebenVect = eigenVects(:,2);
projectedData2 = transpose(activeData) * -nebenVect;
reconstructedData2 = -nebenVect * transpose(projectedData(:,1)) + mean(activeData,2);
plot2DPCA(transpose(activeData),mean(activeData,2),transpose(reconstructedData2),eigenVects,diag(eigenVals),1,1)
plot(projectedData2(:,1),zeros(length(projectedData2(:,1)),2),'+')
title('Task 3b')
meanSquaredError = mean((activeData - reconstructedData2).^2,'all');
%}

%
%   Aufgabe 4
%

activeData = data3D.data;
%activeData = activeData - mean(activeData,2);
[eigenVects3D,eigenVals3D]= ourPca(activeData);
%ourcov = ourCov(activeData);

plot3DPCA(transpose(activeData), transpose(mean(activeData,2)), eigenVects3D, diag(eigenVals3D), 1, 1)

%helperProjection = projectOnHauptVect3D(activeData, eigenVects3D(:,1), eigenVects3D(:,2));
%flip eigenvector sign to get representative data 
projectedData3D = transpose(activeData) * -eigenVects3D;
projectedData3D(:,3) = zeros(1,length(projectedData3D(:,1)));
%plot3(projectedData3D(:,1),projectedData3D(:,2),projectedData3D(:,3),'o');

reconstructedData3D = (-eigenVects3D(:,1:2)) * transpose(projectedData3D(:,1:2)) + mean(activeData,2);
%plot3(reconstructedData3D(1,:),reconstructedData3D(2,:),reconstructedData3D(3,:),'+');
axis equal
title('Task 4')


%
%   Aufgabe 5
%

shapes = dataS.aligned;

%draw all bones 
%{
plot(shapes(:,1,1),shapes(:,2,1));
axis equal
hold on
for i = 2:13
    plot(shapes(:,1,i),shapes(:,2,i));
end
hold off
%}

formattedData = format(shapes);
meanShape = mean(formattedData,1);
formattedData = formattedData - meanShape;
[sortedEigVectsShapes, sortedEigValsShapes] = ourPca(transpose(formattedData));

%use only some eigenvectors
nEigenvectors = 13;

%{
%set b to sqrt(nEigenvalues)
figure
b = 3 * sqrt(sortedEigValsShapes(1:nEigenvectors,1));
temp = b(1,1);
b(:,1) = 0;
b(1,1) = temp;

plotShape(meanShape,sortedEigVectsShapes,b,'b');
axis equal
hold on

b = -b;

plotShape(meanShape,sortedEigVectsShapes,b,'b');
hold off
title('Task 5b')
%}

%100 % of the total variance
temp = sqrt(sortedEigValsShapes(1:nEigenvectors,1));

figure
axis equal
hold on
%plot 5 random example shapes
scale = [0.2, 0.5, 1.0, 2.0, 4.0];
%scale = ones(1,5) * 0.3;
rotation = [30, 60, 90, 120, 0];
%rotation = zeros(1,5);
x_translation = [0, 20, 40, 80, 160];
y_translation = [0,10,20,30,40];
for i = 1:5
    b = (rand(nEigenvectors,1)) .* (temp);
    plotShape(meanShape,sortedEigVectsShapes,b,scale(i),rotation(i),x_translation(i),y_translation(i),'b');
end

hold off
title('Task 5c')



%{
%restrict to 95% of the total variance
sumEigenValues = sum(sortedEigValsShapes);
threshold = 0.95 * sumEigenValues;
A = ones(length(sortedEigValsShapes(:,1)),1) * threshold;
[minValue,closestIndex] = min(abs(cumsum(sortedEigValsShapes)-A));
nEigenvectors = closestIndex;
temp = sqrt(sortedEigValsShapes(1:nEigenvectors,1));

%plot 5 random example shapes
for i = 1:5
b = (rand(nEigenvectors,1)) .* (temp);
plotShape(meanShape,sortedEigVectsShapes,b,'b');
end

hold off
title('Task 5c')
figure
axis equal
hold on
%}

%
%   Util
%

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

% projects Data to an hauptvector
function data = projectOnHauptVect(data, hauptVect)
   for i = 1:length(data)
   data(:,i) = sum(data(:,i).*hauptVect)*hauptVect;
   end
end

function data = projectOnHauptVect3D(data, hauptVect1, hauptVect2)
   for i = 1:length(data)
   item = data(:,i);
   data(:,i) = sum(item.*hauptVect1)*hauptVect1 + sum(item.*hauptVect2)*hauptVect2;
   end
end

