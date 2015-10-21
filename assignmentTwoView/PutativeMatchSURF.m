function [idxs] = PutativeMatchSURF(features1,features2)
%
%  features1 matrix of size no-features x 64 which in the i-th row has 
%            the values of the descriptor of the i-th feature of the first
%            image
%  features2 the feature descriptors for the second image
%
%  idxs      matrix of size no matches x 2 with each row consisting of the
%            feature index for the first image in the first column and the 
%            matching features index in the second image. If a feature does
%            not match no row is provided for the feature


% 
idxs = [];

len1=size(features1,1);
len2=size(features2,1);
K=2;
[D,I] = pdist2(features2,features1,'euclidean','Smallest',K);%K*NoOfFeature1

%but the problem is which should be seen as matched, which doesn't, what's
%the standard: here I choose to compare the most matched two

flag=0.5;%used to judge if match or not
matchRatio=D(1,:)./D(2,:);
ind=find(matchRatio<=flag);
idxs=[ind;I(1,ind)];
idxs=idxs';
end