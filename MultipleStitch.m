function Pano = MultipleStitch( IMAGES, TRANS, fileName )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MultipleStitch 
%   This function stitches multiple images together and outputs the panoramic stitched image
%   with a chain of input images and its corresponding transformations. 
%   Given a chain of images:
%       I1 -> I2 -> I3 -> ... -> Im
%   and its corresponding transformations:
%       T1 transforms I1 to I2
%       T2 transforms I2 to I3 
%       ....
%       Tm-1 transforms Im-1 to Im
%
% INPUT:IMAGES, TRANS, fileName
% OUTPUT:
%   Pano: the final panoramic image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fileName', 'var'),
    fileName = 'pano.jpg';
end

if length(IMAGES) ~= length(TRANS)+1,
    error('Number of images does not match the number of transformations.');
end

outBounds = zeros(2,2);   %去界限
outBounds(1,:) = Inf;
outBounds(2,:) = -Inf;
%选择参照点
refIdx = ceil(median(1:length(IMAGES)));
%计算输出全景图的估计大小
[nrows, ncols, ~] = size(IMAGES{1});
nrows = length(IMAGES) * nrows;
ncols = length(IMAGES) * ncols;

imageToRefTrans = cell(1, length(IMAGES));

%遍历长度，归一化
for idx = 1:length(imageToRefTrans)
    imageToRefTrans{idx} = eye(3);
end

%计算出合适的转换
%左边起
for idx = refIdx-1:-1:1,
    imageToRefTrans{idx} = makeTransformToReferenceFrame(TRANS, idx, refIdx);
    T = imageToRefTrans{idx};
    tmpBounds = findbounds(maketform('affine', T'), [1 1; ncols nrows]);
    outBounds(1,:) = min(outBounds(1,:),tmpBounds(1,:));
    outBounds(2,:) = max(outBounds(2,:),tmpBounds(2,:));
end
%右边起
for idx = refIdx + 1 : length(imageToRefTrans),  
    imageToRefTrans{idx} = makeTransformToReferenceFrame(TRANS, idx, refIdx);
    T = imageToRefTrans{idx};
    T(3, :) = [0, 0, 1]; % Fix rounding errors in the last row.
    tmpBounds = findbounds(maketform('affine', T'), [1 1; ncols nrows]);
    outBounds(1,:) = min(outBounds(1,:),tmpBounds(1,:));
    outBounds(2,:) = max(outBounds(2,:),tmpBounds(2,:));
end

% 拼接Iref图像
XdataLimit = round(outBounds(:,1)');
YdataLimit = round(outBounds(:,2)');
Pano = imtransform( im2double(IMAGES{refIdx}), maketform('affine', eye(3)), 'bilinear', ...
                    'XData', XdataLimit, 'YData', YdataLimit, ...
                    'FillValues', NaN, 'XYScale',1);
             
% 进行转换
%左边起
for idx = refIdx-1:-1:1,
    T = imageToRefTrans{idx};
    Tform = maketform('affine', T');
    AddOn = imtransform(im2double(IMAGES{idx}), Tform, 'bilinear', ...
                        'XData', XdataLimit, 'YData', YdataLimit, ...
                        'FillValues', NaN, 'XYScale',1);
    result_mask = ~isnan(Pano(:,:,1));
    temp_mask = ~isnan(AddOn(:,:,1));
    add_mask = temp_mask & (~result_mask);

    for c = 1 : size(Pano,3),
        cur_im = Pano(:,:,c);
        temp_im = AddOn(:,:,c);
        cur_im(add_mask) = temp_im(add_mask);
        Pano(:,:,c) = cur_im;
    end
end
%右边起
for idx = refIdx + 1 : length(imageToRefTrans),
    T = imageToRefTrans{idx};
    T(3, :) = [0, 0, 1]; % Fix rounding errors in the last row.
    Tform = maketform('affine', T');
    AddOn = imtransform(im2double(IMAGES{idx}), Tform, 'bilinear', ...
                        'XData', XdataLimit, 'YData', YdataLimit, ...
                        'FillValues', NaN, 'XYScale',1);
    result_mask = ~isnan(Pano(:,:,1));
    temp_mask = ~isnan(AddOn(:,:,1));
    add_mask = temp_mask & (~result_mask);

    for c = 1 : size(Pano,3),
        cur_im = Pano(:,:,c);
        temp_im = AddOn(:,:,c);
        cur_im(add_mask) = temp_im(add_mask);
        Pano(:,:,c) = cur_im;
    end
end

%留出未匹配的黑色区域
[I, J] = ind2sub([size(Pano, 1), size(Pano, 2)], find(~isnan(Pano(:, :, 1))));
upper = max(min(I)-1, 1);
lower = min(max(I)+1, size(Pano, 1));
left = max(min(J)-1, 1);
right = min(max(J)+1, size(Pano, 2));
Pano = Pano(upper:lower, left:right,:);

imwrite(Pano, fileName);

end

function T = makeTransformToReferenceFrame(i_To_iPlusOne_Transform, currentFrameIndex, refFrameIndex)
%makeTransformToReferenceFrame
% INPUT:i_To_iPlusOne_Transform, currentFrameIndex, refFrameIndex
%
% OUTPUT:T, frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 YOUR CODE HERE: Calculate T as defined above.                %

if currentFrameIndex<refFrameIndex
    T=eye(3);
    for i=currentFrameIndex:refFrameIndex-1
        T=T*i_To_iPlusOne_Transform{i};
    end
elseif currentFrameIndex>refFrameIndex
    T=eye(3);
    for i=currentFrameIndex-1:refFrameIndex
        T=T*pinv(i_To_iPlusOne_Transform{i});
    end
else
    T=eye(3);
end

%                               END OF YOUR CODE                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
