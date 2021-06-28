clear all

data = load('handdata.mat');
image = data.images{1,1};
%imshow(image);

features = computeFeatures(image);

imagesc(features);
temp = reshape(features(4,:),size(image));
temp = normalize(temp,'range',[0 255]);
imshow(uint8(temp));

function features = computeFeatures(image)
    sz = size(image);
    cols = sz(1)*sz(2);
    features = zeros([8 cols]);
    features(1,:) = reshape(image,[1 cols]);
    
    [FX,FY] = gradient(double(image));
    features(2,:) = reshape(FX,[1 cols]);
    features(3,:) = reshape(FY,[1 cols]);
    
    features(4,:) = sqrt(features(2,:).^2 + features(3,:).^2);
    
    haarlike = computeHaarLike(image);
    features(5,:) = reshape(haarlike(1,:),[1 cols]);
    
    magHaarlike = computeHaarLike(reshape(features(4,:),sz));
    features(6,:) = reshape(magHaarlike(1,:),[1 cols]);
    
    x_coord = (1:sz(1));
    temp = x_coord;
    for i=1 : sz(2)-1
       temp = cat(2,temp,x_coord); 
    end
    features(7,:) = temp;
    
    y_coord = (1:sz(2));
    temp = y_coord;
    for i=1 : sz(1)-1
       temp = cat(2,temp,y_coord); 
    end
    features(8,:) = temp;
    
end
