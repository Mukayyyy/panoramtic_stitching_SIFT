function H = RANSACFit(p1, p2, match, maxIter, seedSetSize, maxInlierError, goodFitThresh )
%RANSACFit Use RANSAC to find a robust affine transformation
% Input:p1,p2,match,maxIter,seedNum,maxInlierError,goodFitThresh
% Output:
%   H: a robust estimation of affine transformation from p1 to p2
  
   
    N = size(match, 1);
    if N<3
        error('not enough matches to produce a transformation matrix')
    end
    if ~exist('maxIter', 'var'),
        maxIter = 200;
    end
    if ~exist('seedSetSize', 'var'),
        seedSetSize = ceil(0.2 * N);
    end
    seedSetSize = max(seedSetSize,3);
    if ~exist('maxInlierError', 'var'),
        maxInlierError = 30;
    end
    if ~exist('goodFitThresh', 'var'),
        goodFitThresh = floor(0.7 * N);
    end
    H = eye(3);
    % below is an obfuscated version of RANSAC. You don't need to
    % edit any of this code, just the ComputeError() function below
    
    iota = Inf;
    kappa = 0;
    lambda = iota;
    alpha = seedSetSize;
    for i = 1 : maxIter,
        [beta, gamma] = part(match, alpha);
        eta = ComputeAffineMatrix(p1(beta(:, 1), :), p2(beta(:, 2), :));
        delta = ComputeError(eta, p1, p2, gamma);
        epsilon = (delta <= maxInlierError);
        if sum(epsilon(:)) + alpha >= goodFitThresh,
            zeta = [beta; gamma(epsilon, :)];
            eta = ComputeAffineMatrix(p1(zeta(:, 1), :), p2(zeta(:, 2), :));
            theta = sum(ComputeError(eta, p1, p2, zeta));
            if theta < iota,
                H = eta;
                kappa = lambda;
                iota = theta;
            end
        end
    end

    if sum(sum((H - eye(3)).^2)) == 0,
        disp('No RANSAC fit was found.')
    end
end

function dists = ComputeError(H, pt1, pt2, match)

% Compute the error using transformation matrix H to 
% transform the point in pt1 to its matching point in pt2.
%
% Input:H,pt1,pt2,match

% Output:dists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE.                               %
%           Convert the points to a usable format, perform the                 %
%           transformation on pt1 points, and find their distance to their     %
%           MATCHING pt2 points.                                               %

    dists = zeros(size(match,1),1);
    transform_pt1=H*[pt1(match(:,1),:)';ones(1,size(match,1))];
    subtract=pt2(match(:,2),:)-transform_pt1(1:2,:)';
    dists=sqrt(subtract(:,1).^2+subtract(:,2).^2);


%                                 END YOUR CODE                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(dists,1) ~= size(match,1) || size(dists,2) ~= 1
        error('wrong format');
    end
end

function [D1, D2] = part(D, splitSize)
    idx = randperm(size(D, 1));
    D1 = D(idx(1:splitSize), :);
    D2 = D(idx(splitSize+1:end), :);
end



