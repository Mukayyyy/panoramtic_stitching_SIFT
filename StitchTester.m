clc; close all; clc;
imgList = dir('./data/grass*.jpg');  %����һϵ��ͼƬ
saveFileName = 'yosemite.jpg';       %��ƴ�Ӻ��ͼ��
addpath('KeypointDetect');           %��·��
IMAGES = cell(1, length(imgList));

for i = 1 : length(imgList)
    IMAGES{i} = imread(['./data/' imgList(i).name]);
    if max(size(IMAGES{i})) > 1000 || length(imgList) > 10
        IMAGES{i} = imresize(IMAGES{i}, 0.6);    %��ֹͼ��̫���ʱ̫�� ���¶���ͼ���С
    end
end
disp('Images loaded. Beginning feature detection...');

%���ú���SIFTDescriptor���SIF������
DESCRIPTOR = cell(1, length(imgList));
POINT_IN_IMG = cell(1, length(imgList));
for i = 1 : length(imgList)
    [feature, ~, imp] = detect_features(IMAGES{i});
    POINT_IN_IMG{i} = feature(:, 1:2);
    pointInPyramid = feature(:, 8:9);
    DESCRIPTOR{i} = SIFTDescriptor(imp, pointInPyramid, feature(:,3));
end

%����RANSACFit�����任affine����
TRANSFORM = cell(1, length(imgList)-1);
for i = 1 : (length(imgList)-1)
    disp(['fitting transformation from ' num2str(i) ' to ' num2str(i+1)])
    M = SIFTSimpleMatcher(DESCRIPTOR{i}, DESCRIPTOR{i+1}, 0.7);
    TRANSFORM{i} = RANSACFit(POINT_IN_IMG{i}, POINT_IN_IMG{i+1}, M);
end

%����ȫ��ͼ�������
disp('Stitching images...')
MultipleStitch(IMAGES, TRANSFORM, saveFileName);
disp(['The completed file has been saved as ' saveFileName]);
imshow(imread(saveFileName));