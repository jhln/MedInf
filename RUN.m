
clear all

%data load
data2D = load('daten.mat');
data3D = load('daten3d.mat');
dataS = load('shapes.mat');



%x = data2D.data1(1,:);
%y = data2D.data1(2,:);

%plot(x,y,'*');
%plot(data2D.data1(1,:), data2D.data1(2,:))
%axis equal


%
%   Aufgabe 5
%

dataShapes = dataS.aligned;
meanShape = mean(dataShapes,3);
[sortedEigVecS, sortedEigValS] = pca(dataShapes);
%sortedEigVecS(:,2)
chosenEigenV = [sortedEigVecS(:,1), sortedEigVecS(:,2),sortedEigVecS(:,3) ];
dummyReconstructionS = generateShapes([1,1,1], chosenEigenV )

% noch zu schreiben, plottet die reconstructed Shapes anhand der varianc
% und abweichung 
%%%%plotShape(dataShapes, dmeanShape, dummyReconstructionS, sortedEigVecS, sortedEigValS, 1, 1)



%
%   Aufgabe 4
%

data3 = data3D.data;
mean(data3,2);
%dataMean3 = data3 - mean(data3,2);
[sortedEigVec3,sortedEigVal3]= pca(data3); % meint ihr das ist realistisch?
%mean(data3,2)
%diag(sortedEigVal3)
%%%%plot3DPCA(transpose(data3), transpose(mean(data3,2)), sortedEigVec3, diag(sortedEigVal3), 1, 1)

% error elipsoid fehlt noch

%und und ob das jetzt stimmt ... bin mir nicht sicher, wie man auf eine
%fläche projeziert, also mit 2 Hauptvectoren - ein offener versuch unter projectOn2HauptVect
%Cloud3D*[sortedEigVec3(:,1),sortedEigVec3(:,2)]
projected3 = projectOn2HauptVect(data3, sortedEigVec3(:,1), sortedEigVec3(:,2));
%%%%plot3DPCA(transpose(projected3), transpose(mean(projected3,2)), sortedEigVec3, diag(sortedEigVal3), 1, 1)



%
%   Aufgabe 3
%

data2 = data2D.data2;
dataMean2 = data2 - mean(data2,2);

[sortedEigVec2, sortedEigVal2] = pca(data2);
hauptVect2 = sortedEigVec2(:,1);
%plot(d2(1,:),d2(2,:))
dummyReconstruction2 = projectOnHauptVect(data2, hauptVect2);

% nicht sicher wegen den parametern - die daten müssten meiner meinung nach
% transponiert wernden, but not sure
%%%%plot2DPCA(transpose(data2),mean(d2,2),transpose(projectedData2),eigenVects2,diag(eigenVals2),1,1)
%bin mir auch nicht sicher, was der richtige mean wert ist
%%%%plot2DPCA(transpose(data2),dataMean2,transpose(dummyReconstruction2),sortedEigVec2,diag(sortedEigVal2),1,1)
%%%%plot2DPCA(data2, dataMean2, dummyReconstruction2, sortedEigVec2, sortedEigVal2, 1, 1)

%nochmal für den 1. Nebenvektor - wobei es ja eigentlich nur 2 überhaupt
%gibt !?
nebenVect2 = sortedEigVec2(:,2);
%plot(d2(1,:),d2(2,:))
projectedData2 = projectOnHauptVect(data2, nebenVect2);

%%%%plot2DPCA(transpose(d2),mean(d2,2),transpose(projectedData2),eigenVects2,diag(eigenVals2),1,1)


%
%   Aufgabe 2
%

d1 = data2D.data1;

[eigenVects1, eigenVals1] = pca(d1);
hauptVect1 = eigenVects1(:,1);
%plot(d2(1,:),d2(2,:))

% nicht sicher wegen den parametern - weil hier keine reconstruction
% stattfindet
%%%%plot2DPCA(transpose(d1),mean(d1,2),transpose(d1),eigenVects1,diag(eigenVals1),1,0)


%plotShape(dataShapes, dmeanShape, dummyReconstructionS, sortedEigVecS, sortedEigValS, 1, 1)



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



function data = projectOn2HauptVect (data, hauptVect1, hauptVect2)
   for i = 1:length(data)
   %data(:,i)
   item = data(:,i);
   data(:,i) = sum(item.*hauptVect1)*hauptVect1 + sum(item.*hauptVect2)*hauptVect2;
   %data(:,i)
   end
   return
end


% projects Data to an hauptvector
function data = projectOnHauptVect (data, hauptVect)
   for i = 1:length(data)
   %data(:,i)
   data(:,i) = sum(data(:,i).*hauptVect)*hauptVect;
   %data(:,i)
   end
   return
end

%computes pca
function [geordneteEigenVects, geordneteEigenVals] = pca(data)
    %[v,d] = eig(data)
    
    data = data - mean(data,2);
    %nexttile
    %plot(data(1,:),data(2,:))
    %[a,b] = eig(ourCov2(data))
    [geordneteEigenVals,i] = sort(eig(ourCov2(data)),'descend');
    [vect,val] = eig(ourCov2(data));
    geordneteEigenVects = vect(:,i);
    return
end

%computes covariance
function covMatrix = ourCov2(dataMatrix)

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
    %data = load('/home/johannes/uni/21ss/MedInf/daten.mat');
    %m = cov(data.data1(:,1),data.data1(:,2));
    %disp(m);
    
    weight = 1/(length(column1) - 1 );
    
    mean1 = mean(column1);
    mean2 = mean(column2);
    
    result = weight * (sum((column1 - mean1) .* (column2 - mean2)));
    return
end

%{
[eigenValues,i] = sort(eig(ourCov2(d)),'descend')
[vect,val] = eig(ourCov2(d))
geordneteEigenVects = vect(:,i)
geordneteEigenVects(:,1)
%}
%hauptVect = geordneteEigenVects(:,1)

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


%{
fn = fieldnames(data2D)
c = numel(fn)

tiledlayout(c,1)

for i = 1:c
    d = data2D.(fn{i});
    d = d - mean(d,2)

    nexttile
    plot(d(1,:), d(2,:),'x')
    axis equal

    title('data')
    cov(d(1,:),d(2,:))
    ourCov2(d)
end
%}

%task 1a with own covFunction
%ownMat2D = ourCov2(data2D.data1)
%task 1a comparison with predefined cov
%m2D = cov(data2D.data1(1,:),data2D.data1(2,:))
