function match = SIFTSimpleMatcher(descriptor1, descriptor2, thresh)
% SIFTSimpleMatcher 
%   Match one set of SIFT descriptors (descriptor1) to another set of
%   descriptors (decriptor2). Each descriptor from descriptor1 can at
%   most be matched to one member of descriptor2, but descriptors from
%   descriptor2 can be matched more than once.
% INPUT:descriptor1,descriptor2,thresh
% OUTPUT: Match(N * 2 matrix, each row is a match)

    if ~exist('thresh', 'var'),
        thresh = 0.7;
    end

    match = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE:                               %

[N1,~]=size(descriptor1);
[N2,~]=size(descriptor2);
for i=1:N1
    distance=[];
    for j=1:N2
        subtract=descriptor1(i,:)-descriptor2(j,:);
        distance=[distance,norm(subtract)];
    end
        sort_distance=sort(distance);
        if(sort_distance(1)<0.7*sort_distance(2))
            j=find(distance==sort_distance(1));
            match=[match;[i,j]];
        end
end
%                                 END YOUR CODE                                %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
