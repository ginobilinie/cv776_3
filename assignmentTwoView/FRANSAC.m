function [fMat,inliers]=FRANSAC(matchedPoints1, matchedPoints2,probSol,threshold)
%
%  matchedPoints1, matchedPoints2 matrix of size no_matches x 2 containing
%                                 the x, y coordinates of all matched 
%                                 features in image 1,2 respectively
%
%  probSol  probability of having seen a good solution before stopping. 
%           Typically,this should be 0.95 to 0.99
%  
%  threshold inlier threshold typical ranges between 1 und 4 pixels
%
%

% 

fMat=eye(3);
inliers=0;

%to choose random samples
n=8;
len1=length(matchedPoints1);
len2=length(matchedPoints2);
if len1~=len2
    fprintf('two matched points have different length\n');
end
chosenInd=randperm(len1,n);

%# of iterations
e=0.5;
N=round(log(1-probSol)/log(1-(1-e)^8));

num1=size(matchedPoints1,1);
num2=size(matchedPoints2,1);
if num1~=num2
    fprintf('two matched points have different length\n');
end
matchedPoints1=[matchedPoints1,ones(num1,1)];
matchedPoints2=[matchedPoints2,ones(num1,1)];

%main steps of RANSAC
optimalInliers=[];
for i=1:N
    %pick 4 pairs randomly
    randInds=randi(num1,8,1);
    randPoints1=matchedPoints1(randInds,:);
    randPoints2=matchedPoints2(randInds,:);
    
    % Hartley normalization
    normMat1=HartleyNorm(randPoints1);
    normMat2=HartleyNorm(randPoints2);
    randPoints1=(normMat1*randPoints1');
    randPoints2=(normMat2*randPoints2');
    
    %computer foundamental matrix
    fMat=genFoundamentalMat(randPoints1,randPoints2,normMat1,normMat2);
   
    currInliers=[];
    for k=1:num1
        p1=matchedPoints1(k,:)';
        p2=matchedPoints2(k,:)';
        tmp=p2' * fMat * p1;
        % threshold check
        if(abs(tmp)<=threshold)
            currInliers(end+1)=k;
        end      
    end
    
    %keep the largest set of inliers
    if(size(currInliers,2)>size(optimalInliers,2))
        optimalInliers=currInliers;
    end
end

%recomputer fMat
inliers=optimalInliers';
randPoints1=matchedPoints1(inliers,:);
randPoints2=matchedPoints2(inliers,:);
normMat1=HartleyNorm(randPoints1);
normMat2=HartleyNorm(randPoints2);
randPoints1=(normMat1*randPoints1');
randPoints2=(normMat2*randPoints2');
fMat=genFoundamentalMat(randPoints1,randPoints2,normMat1,normMat2);
return

%Hartley normalization
function hMat=HartleyNorm( points )
num1=size(points,1);
points=points';
center=mean(points,2);
dist=mean(sqrt(sum((points - repmat(center,1,num1)).^2,1)));
hMat=[sqrt(2)/dist, 0, -sqrt(2)/dist*center(1); 0, sqrt(2)/dist, -sqrt(2)/dist*center(2); 0, 0, 1];
return

%generate foundamental mat
function fMat=genFoundamentalMat(points1,points2,norm1,norm2)
A=[repmat(points2(1,:)',1,3) .* points1', repmat(points2(2,:)',1,3) .* points1', points1'];
[U,S,V]=svd(A);
mat=reshape(V(:,end),3,3)';
% constrain
[u,s,v]=svd(mat);
matPrime=u*diag([s(1) s(5) 0])*(v');
fMat= norm2' * matPrime* norm1;
return










