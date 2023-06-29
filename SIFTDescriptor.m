function descriptors = SIFTDescriptor(pyramid, keyPtLoc, keyPtScale)
% SIFTDescriptor Build SIFT descriptors from image at detected key points'
% location with detected key points' scale and angle
% ����SIF�����Ӻ��� ���ڼ��ؼ���λ�õı����ͽǶ�
% INPUT:pyramid: Image pyramid;  keyPtLoc; keyPtScale; 
% OUTPUT:descriptors: N * 128 matrix, each row is a feature descriptor

    grad_theta = cell(length(pyramid),1);
    grad_mag = cell(length(pyramid),1);
    
    for scale = 1:length(pyramid)
        currentImage = pyramid{scale};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                YOUR CODE HERE                                %

        % x�����ݶ�
        img_dx = zeros(size(currentImage));
        img_dx=filter2([-1 0 1],currentImage);   %ʹ������ͼ���С����
        % y�����ݶ�
        img_dy = zeros(size(currentImage));
        img_dy=filter2([-1;0;1],currentImage);
        
        %����xy�����ϵ��ݶ�ֵ
        grad_mag{scale} = zeros(size(currentImage));
        grad_theta{scale} = zeros(size(currentImage));
        grad_mag{scale}=sqrt(img_dx.^2+img_dy.^2);
        grad_theta{scale}=atan2(img_dy,img_dx);   
        
%                               END OF YOUR CODE                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        grad_theta{scale} = mod(grad_theta{scale}, 2*pi);  %���㷽��
        
    end
    
    %���ýǶ��ݶȵĲ���ֵ
    num_angles = 8;
    num_histograms = 4; 
    pixelsPerHistogram = 4;
    
    patch_size = num_histograms * pixelsPerHistogram;  %��ȡ���ؼ��㸽���Ŀռ��С
    
    N = size(keyPtLoc, 1); 
    descriptors = zeros(N, num_histograms * num_histograms * num_angles);
            
    % ��������ÿһ���ؼ���
    for i = 1 : N
        scale = round(keyPtScale(i));    
        %Ѱ�ҹؼ��������
        %ȷ��DOG�ؼ��������
        xAtScale = keyPtLoc(i, 1);   
        yAtScale = keyPtLoc(i, 2);
        x_lo = round(xAtScale - patch_size / 2);
        x_hi = x_lo+patch_size-1;
        y_lo = round(yAtScale - patch_size / 2);
        y_hi = y_lo+patch_size-1;
                
        magnitudes = grad_mag{scale};
        thetas = grad_theta{scale};
        try    
            patch_mag = zeros(patch_size,patch_size);   %��ȡɸѡ�ؼ���
            patch_theta = zeros(patch_size,patch_size);
            patch_mag = magnitudes(y_lo:y_hi,x_lo:x_hi);
            patch_theta = thetas(y_lo:y_hi,x_lo:x_hi);
        catch err
            continue;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              YOUR CODE HERE:                                 %

        % ���������ݶȷ���
        patch_angle_offset = ComputeDominantDirection(patch_mag, patch_theta);
        patch_theta = patch_theta-patch_angle_offset;

        patch_theta = mod(patch_theta, 2*pi);  %�ض�λ
        patch_mag = patch_mag .* fspecial('gaussian', patch_size, patch_size / 2);  %ʹ�ø�˹��������
                 
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
       
        descriptors(i, :) = feature;   %�������������������
    end
    descriptors = NormalizeDescriptors(descriptors);    %��һ��������
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

    %�����ݶȵı仯��
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
    %����ֱ��ͼ
    [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles);
        %�������ֵ��ȷ��ֱ��ͼ�ķ���
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
% INPUT��descriptors(N x 128 matrix where each row is a SIFT descriptor)
% OUTPUT��descriptors(N x 128 matrix containing a normalized version of the
% input)
    
    %��һ��
    lengths = sqrt(sum(descriptors.^2, 2));   
    nonZeroIndices = find(lengths);
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

    descriptors(descriptors > 0.2) = 0.2;
    lengths = sqrt(sum(descriptors.^2, 2));
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

end