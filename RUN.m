
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

activeData = activeData - mean(activeData,2);

%plot(activeData(1,:), activeData(2,:),'x')
%axis equal
%title('Aufgabe1')
 


%
%   Aufgabe 2
%

[eigenVects, eigenVals] = ourPca(activeData);
hauptVect = eigenVects(:,1);
%plot2DPCA(transpose(activeData),mean(activeData,2),0,eigenVects,diag(eigenVals),1,0)
%title('Aufgabe2')



%
%   Aufgabe 3
%
%shift data around mean done in task1. 

[eigenVects, eigenVals] = ourPca(activeData);
hauptVect = eigenVects(:,1);
helperProjection = projectOnHauptVect(activeData, hauptVect);

%andere art zu projezieren siehe Task 4, auch hier scheinen die daten
%um die x-achse gespiegelt, vielleicht wegen dem transpose?
projectedData = transpose(activeData) * eigenVects;
projectedData(:,2) = zeros(1,length(projectedData(:,2)));
reconstructedData = (hauptVect) * transpose(projectedData(:,1)) + mean(activeData,2);
%additionally ploting the correct helperProjection (o) and the
%reconstruction (+)
plot2DPCA(transpose(activeData),mean(activeData,2),(projectedData),eigenVects,diag(eigenVals),1,1);
plot(helperProjection(1,:),helperProjection(2,:),'o')
plot(reconstructedData(1,:),reconstructedData(2,:),'+')
title('Aufgabe3.1')

nebenVect = eigenVects(:,2);
projectedData2 = projectOnHauptVect(activeData, nebenVect);
%plot2DPCA(transpose(activeData),mean(activeData,2),transpose(projectedData2),eigenVects,diag(eigenVals),1,1)
%title('Aufgabe3.2')



%
%   Aufgabe 4
%

activeData = data3D.data;
activeData = activeData - mean(activeData,2);
[eigenVects3D,eigenVals3D]= ourPca(activeData);

plot3DPCA(transpose(activeData), transpose(mean(activeData,2)), eigenVects3D, diag(eigenVals3D), 1, 1)
% error elipsoid fehlt noch

%projectedData3D = projectOnHauptVect3D(activeData, eigenVects3D(:,1), eigenVects3D(:,2));

%andere art zu projezieren, nicht ganz sicher obs stimmt, z-komponente sollte
%vernachlässigbar sein, daten wirken aber um x gespiegelt, aber die rekonstruktion scheint richtig.
projectedData3D = transpose(activeData) * eigenVects3D;
projectedData3D(:,3) = zeros(1,length(projectedData3D(:,1)));
plot3(projectedData3D(:,1),projectedData3D(:,2),projectedData3D(:,3),'o');

reconstructedData3D = (eigenVects3D(:,1:2)) * transpose(projectedData3D(:,1:2)) + mean(activeData,2);
plot3(reconstructedData3D(1,:),reconstructedData3D(2,:),reconstructedData3D(3,:),'+');
title('Aufgabe 4')

return


%
%   Aufgabe 5
%

dataShapes = dataS.aligned;
meanShape = mean(dataShapes,3);
[sortedEigVecS, sortedEigValS] = ourPca(dataShapes);
%sortedEigVecS(:,2)
chosenEigenV = [sortedEigVecS(:,1), sortedEigVecS(:,2),sortedEigVecS(:,3) ];
dummyReconstructionS = generateShapes([1,1,1], chosenEigenV )

% noch zu schreiben, plottet die reconstructed Shapes anhand der varianc
% und abweichung 
%%%%plotShape(dataShapes, dmeanShape, dummyReconstructionS, sortedEigVecS, sortedEigValS, 1, 1)



%
%   Util
%

generateShapes([1, 0.5, 0.1], [ 1,2,3;1,2,3;1,2,3;1,2,3;]);

function generatedShape = generateShapes(paraVec, eigenVects)
    l = length(eigenVects);
    generatedShape = transpose([zeros(1,l)]);
    for i = 1:length(paraVec)
        %i
        %p = paraVec(i)
        %eigenVects(:,i)
        %paraVec(i) * eigenVects(:,i)
        generatedShape = generatedShape + (paraVec(i) * eigenVects(:,i));
    end
    return
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


%d(:,1)
%{
for i = 1:length(d)
   d(:,i)
   d(:,i) = sum(d(:,i).*hauptVect)*hauptVect;
   d(:,i)
end
%}

%hauptVect
%sum(d(:,1).*hauptVect).*hauptVect


%@(x) sum(d(:,x).*hauptVect)*hauptVect
%C=(sum(A.*B)/(norm(B)^2))*B;
%sum([1 0 0] .* [1 0 0])


