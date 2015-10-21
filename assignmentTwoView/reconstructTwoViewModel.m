% Authors: Johannes L. Schönberger, Jan-Michael Frahm <{jsch,jmf} at cs.unc.edu>

function [projMatrix1, projMatrix2, points3D] = reconstructTwoViewModel(image1file, image2file, visualize)

% load image data (including features, exif, gray scale image)
load(image1file);
image1name = imagename;
image1=image_data;
image1_exif = image_exif;
points1 = points;
features1 = features;


load(image2file);
image2name = imagename;
image2=image_data;
image2_exif = image_exif;
points2 = points;           % feature locations
features2 = features;       % feature descriptors

% get focal length from exif
K1 = getCalibrationFromExif(image1_exif);
K2 = getCalibrationFromExif(image2_exif);

if visualize
    % randomly choose 200 points for visualization
    index1 = randi(size(points1,1),1,200);
    index2 = randi(size(points2,1),1,200);
    
    figure;
    
    subplot(1, 2, 1);
    imshow(image1);
    hold on;
    plot(points1(index1,1),points1(index1,2),'*')
    title('200 random keypoints for image 1');
    
    subplot(1, 2, 2);
    imshow(image2);
    hold on;
    plot(points2(index2,1),points2(index2,2),'*')
    title('200 random keypoints for image 2');
end

indexPairs = matchFeaturesSURF(features1, features2);

matchedPoints1 = points1(indexPairs(:,1),:);
matchedPoints2 = points2(indexPairs(:,2),:);

if visualize
    figure;
    VisualizeMatchedFeatures(image1, image2, matchedPoints1, matchedPoints2);
end


[F, inliers] = FRANSAC(matchedPoints1, matchedPoints2,0.99,2);

inlierPoints1 = matchedPoints1(inliers,:);
inlierPoints2 = matchedPoints2(inliers,:);

% Compose essential matrix using camera calibration
E = K2' * F * K1;

% Normalize image coordinages
normPoints1 = inv(K1) * [inlierPoints1 ones(size(inlierPoints1, 1), 1)]';
normPoints1 = normPoints1(1:2,:)';
normPoints2 = inv(K2) * [inlierPoints2 ones(size(inlierPoints2, 1), 1)]';
normPoints2 = normPoints2(1:2,:)';

[projMatrix1, projMatrix2, points3D] = ...
    relativePoseFromEssentialMatrix(E, normPoints1, normPoints2);


% Filter unstable 3D points
angles = calcTriangulationAngles(projMatrix1, projMatrix2, points3D);
minAngle = 2 * pi / 180;
points3D = points3D(angles > minAngle,:);
depths1 = calcDepths(projMatrix1, points3D);
depths2 = calcDepths(projMatrix2, points3D);
points3D = points3D(depths1 > 0 & depths2 > 0,:);

if visualize
    figure;
    plot3(points3D(:,1), points3D(:,2), points3D(:,3), 'b.');
    hold on;
    
    invProjMatrix2 = inv([projMatrix2; 0 0 0 1]);
    
    cameraModel = [-1 -1 0 1
                   -1 1 0 1
                   1 1 0 1
                   1 -1 0 1
                   -1 -1 0 1
                   nan nan 0 1
                   0 0 -1 1
                   -1 -1 0 1
                   nan nan 0 1
                   0 0 -1 1
                   -1 1 0 1
                   nan nan 0 1
                   0 0 -1 1
                   1 1 0 1
                   nan nan 0 1
                   0 0 -1 1
                   1 -1 0 1];
               
    cameraModel(:,1:3) = 0.2 * cameraModel(:,1:3);
               
    cameraModel1 = cameraModel;
    cameraModel2 = invProjMatrix2(1:3,:) * cameraModel';
    cameraModel2 = cameraModel2';
    
    plot3(cameraModel1(:,1), cameraModel1(:,2), cameraModel1(:,3), 'r-');
    plot3(cameraModel2(:,1), cameraModel2(:,2), cameraModel2(:,3), 'g-');
    
    grid on;
    axis equal;
    xlim([-5 5]);
    ylim([-5 5]);
    zlim([-1 11]);
    
    legend('Point cloud', 'Camera 1', 'Camera 2');
end

end

function [K]=getCalibrationFromExif(ex)

        % call external function for getting K from exif
        K=Exif2Calibration(ex);
 
end

function [R1, R2, t1, t2] = decomposeEssentialMatrix(E)

[U, ~, V] = svd(E);

if det(U) < 0
    U = -U;
end
if det(V) < 0
    V = -V;
end

W = [0 -1 0;
     1 0 0;
     0 0 1];

R1 = U * W * V';
R2 = U * W' * V';

t1 = U(:,3);
t2 = -U(:,3);

end

function points3D = triangulatePoints(projMatrix1, projMatrix2, points2D1, points2D2)

numPoints = size(points2D1, 1);

points3D = zeros(numPoints, 3);

for i = 1 : numPoints
    A = zeros(4, 4);
    A(1,1) = projMatrix1(1,1) - points2D1(i,1) * projMatrix1(3,1);
    A(1,2) = projMatrix1(1,2) - points2D1(i,1) * projMatrix1(3,2);
    A(1,3) = projMatrix1(1,3) - points2D1(i,1) * projMatrix1(3,3);
    A(1,4) = projMatrix1(1,4) - points2D1(i,1) * projMatrix1(3,4);
    A(2,1) = projMatrix1(2,1) - points2D1(i,2) * projMatrix1(3,1);
    A(2,2) = projMatrix1(2,2) - points2D1(i,2) * projMatrix1(3,2);
    A(2,3) = projMatrix1(2,3) - points2D1(i,2) * projMatrix1(3,3);
    A(2,4) = projMatrix1(2,4) - points2D1(i,2) * projMatrix1(3,4);
    A(3,1) = projMatrix2(1,1) - points2D2(i,1) * projMatrix2(3,1);
    A(3,2) = projMatrix2(1,2) - points2D2(i,1) * projMatrix2(3,2);
    A(3,3) = projMatrix2(1,3) - points2D2(i,1) * projMatrix2(3,3);
    A(3,4) = projMatrix2(1,4) - points2D2(i,1) * projMatrix2(3,4);
    A(4,1) = projMatrix2(2,1) - points2D2(i,2) * projMatrix2(3,1);
    A(4,2) = projMatrix2(2,2) - points2D2(i,2) * projMatrix2(3,2);
    A(4,3) = projMatrix2(2,3) - points2D2(i,2) * projMatrix2(3,3);
    A(4,4) = projMatrix2(2,4) - points2D2(i,2) * projMatrix2(3,4);

    [~, ~, V] = svd(A);

    points3D(i,1) = V(1,4) / V(4,4);
    points3D(i,2) = V(2,4) / V(4,4);
    points3D(i,3) = V(3,4) / V(4,4);
end

end

function depths = calcDepths(projMatrix, points3D)

points2D = projMatrix * [points3D ones(size(points3D, 1), 1)]';

depths = norm(projMatrix(1:3,3)) * points2D(3,:);

end

function [projMatrix1, projMatrix2, points3D] = ...
    relativePoseFromEssentialMatrix(E, points2D1, points2D2)

[R1, R2, t1, t2] = decomposeEssentialMatrix(E);

projMatrix1 = eye(3, 4);
projMatrices = {[R1 t1], [R1 t2], [R2 t2] [R2 t1]};

bestNumCheralityPasses = 0;
bestProjMatrix2 = nan;
bestPoints3D = nan;

for i = 1 : 4
    projMatrix2 = projMatrices{i};
    
    points3D = triangulatePoints(projMatrix1, projMatrix2, points2D1, points2D2);
    depths = calcDepths(projMatrix2, points3D);
    
    numCheralityPasses = sum(depths > 0);
    
    if numCheralityPasses > bestNumCheralityPasses
        bestNumCheralityPasses = numCheralityPasses;
        bestProjMatrix2 = projMatrix2;
        bestPoints3D = points3D;
    end
end

projMatrix2 = bestProjMatrix2;
points3D = bestPoints3D;

end

function angles = calcTriangulationAngles(projMatrix1, projMatrix2, points3D)

invProjMatrix1 = inv([projMatrix1; 0 0 0 1]);
invProjMatrix2 = inv([projMatrix2; 0 0 0 1]);

baseline = norm(invProjMatrix1(1:3,4) - invProjMatrix2(1:3,4));

numPoints3D = size(points3D, 1);

dists1 = points3D - repmat(invProjMatrix1(1:3,4)', [numPoints3D 1]);
dists1 = sqrt(sum(dists1.^2, 2));
dists2 = points3D - repmat(invProjMatrix2(1:3,4)', [numPoints3D 1]);
dists2 = sqrt(sum(dists2.^2, 2));

angles = abs(acos((dists1.^2 + dists1.^2 - baseline.^2) ./ (2 .* dists1 .* dists2)));

end

function idxs = matchFeaturesSURF(features1, features2)

    idxs = PutativeMatchSURF(features1,features2);

end


function []= VisualizeMatchedFeatures(image1, image2, matchedPoints1, matchedPoints2)
    
    % visualize features
    images = [image1,image2];
    figure, imshow(images), hold on
    index = randi(size(matchedPoints1,1),1,200);

    % plot features
    plot(matchedPoints1(index,1),matchedPoints1(index,2),'*')
    plot(matchedPoints2(index,1)+size(image1,1),matchedPoints2(index,2),'*')

    % plot lines
    for i = 1 : size(index,2)
       plot([matchedPoints1(index(i),1),matchedPoints2(index(i),1)+size(image1,2)],...
           [matchedPoints1(index(i),2),matchedPoints2(index(i),2)]);
    end


end
