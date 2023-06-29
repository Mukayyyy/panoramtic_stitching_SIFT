function descriptors = SIFTDescriptor(pyramid, keyPtLoc, keyPtScale)
% SIFTDescriptor Build SIFT descriptors from image at detected key points'
% location with detected key points' scale and angle
% 构建SIF描述子函数 用于检测关键点位置的比例和角度
% INPUT:pyramid: Image pyramid;  keyPtLoc; keyPtScale; 
% OUTPUT:descriptors: N * 128 matrix, each row is a feature descriptor

    grad_theta = cell(length(pyramid),1);
    grad_mag = cell(length(pyramid),1);
    
    for scale = 1:length(pyramid)
        currentImage = pyramid{scale};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE                                %

        % x方向梯度
        img_dx = zeros(size(currentImage));
        img_dx=filter2([-1 0 1],currentImage);   %使处理后的图像大小不变
        % y方向梯度
        img_dy = zeros(size(currentImage));
        img_dy=filter2([-1;0;1],currentImage);
        
        %计算xy方向上的梯度值
        grad_mag{scale} = zeros(size(currentImage));
        grad_theta{scale} = zeros(size(currentImage));
        grad_mag{scale}=sqrt(img_dx.^2+img_dy.^2);
        grad_theta{scale}=atan2(img_dy,img_dx);   
        
%                               END OF YOUR CODE                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        grad_theta{scale} = mod(grad_theta{scale}, 2*pi);  %计算方便
        
    end
    
    %设置角度梯度的参数值
    num_angles = 8;
    num_histograms = 4; 
    pixelsPerHistogram = 4;
    
    patch_size = num_histograms * pixelsPerHistogram;  %提取各关键点附近的空间大小
    
    N = size(keyPtLoc, 1); 
    descriptors = zeros(N, num_histograms * num_histograms * num_angles);
            
    % 迭代遍历每一个关键点
    for i = 1 : N
        scale = round(keyPtScale(i));    
        %寻找关键点的象限
        %确定DOG关键点的中心
        xAtScale = keyPtLoc(i, 1);   
        yAtScale = keyPtLoc(i, 2);
        x_lo = round(xAtScale - patch_size / 2);
        x_hi = x_lo+patch_size-1;
        y_lo = round(yAtScale - patch_size / 2);
        y_hi = y_lo+patch_size-1;
                
        magnitudes = grad_mag{scale};
        thetas = grad_theta{scale};
        try    
            patch_mag = zeros(patch_size,patch_size);   %提取筛选关键点
            patch_theta = zeros(patch_size,patch_size);
            patch_mag = magnitudes(y_lo:y_hi,x_lo:x_hi);
            patch_theta = thetas(y_lo:y_hi,x_lo:x_hi);
        catch err
            continue;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              YOUR CODE HERE:                                 %

        % 计算主导梯度方向
        patch_angle_offset = ComputeDominantDirection(patch_mag, patch_theta);
        patch_theta = patch_theta-patch_angle_offset;

        patch_theta = mod(patch_theta, 2*pi);  %重定位
        patch_mag = patch_mag .* fspecial('gaussian', patch_size, patch_size / 2);  %使用高斯函数测量
                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             YOUR CODE HERE:                                  %
%       Compute the gradient histograms and concatenate them in the          
%       feature variable to form a size 1x128 SIFT descriptor for this keypoint

        feature = [];
        subdivided_patch_theta=zeros(pixelsPerHistogram,pixelsPerHistogram);
        subdivided_patch_mag=zeros(pixelsPerHistogram,pixelsPerHistogram);
        for y=1:4:13
            for x=1:4:13
                subdivided_patch_theta=patch_theta(y:y+3,x:x+3);
                subdivided_patch_mag=patch_mag(y:y+3,x:x+3);
                [histogram,angles]=ComputeGradientHistogram(num_angles,subdivided_patch_mag,subdivided_patch_theta);
                feature=[feature,histogram];
            end
        end
        
%                                 END YOUR CODE                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        descriptors(i, :) = feature;   %给描述子添加特征向量
    end
    descriptors = NormalizeDescriptors(descriptors);    %归一化描述子
end



function [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles)
% Compute a gradient histogram using gradient magnitudes and directions.
% Each point is assigned to one of num_bins depending on its gradient
% direction; the gradient magnitude of that point is added to its bin.
%
% INPUT: num_bins, gradient_magnitudes, gradient angles:                                     
% OUTPUT: histogram, angles

    angle_step = 2 * pi / num_bins;
    angles = 0 : angle_step : (2*pi-angle_step);

    %测量梯度的变化率
    histogram = zeros(1, num_bins);
    [rows,cols]=size(gradient_magnitudes);
    for m=1:rows
        for n=1:cols
            for angle=0:angle_step:(2*pi-angle_step)
                if(angle<=gradient_angles(m,n)&&gradient_angles(m,n)<angle+angle_step)
                    histogram(round(angle/angle_step)+1)=histogram(round(angle/angle_step)+1)+gradient_magnitudes(m,n);
                    break;
                end
            end
        end
    end
end



function direction = ComputeDominantDirection(gradient_magnitudes, gradient_angles)
% Computes the dominant gradient direction for the region around a keypoint
% given the scale of the keypoint and the gradient magnitudes and gradient
% angles of the pixels in the region surrounding the keypoint.
%
% INPUT:gradient_magnitudes, gradient_angles, scale
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE:                               %

    num_bins = 36;
    %计算直方图
    [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles);
        %计算最大值，确定直方图的方向
        peak=max(histogram);
        loc=find(histogram==peak);
        direction =2*pi*loc(1)/num_bins;

%                                 END YOUR CODE                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



function descriptors = NormalizeDescriptors(descriptors)
% Normalizes SIFT descriptors so they're unit vectors. You don't need to
% edit this function.
%
% INPUT：descriptors(N x 128 matrix where each row is a SIFT descriptor)
% OUTPUT：descriptors(N x 128 matrix containing a normalized version of the
% input)
    
    %归一化
    lengths = sqrt(sum(descriptors.^2, 2));   
    nonZeroIndices = find(lengths);
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

    descriptors(descriptors > 0.2) = 0.2;
    lengths = sqrt(sum(descriptors.^2, 2));
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

end