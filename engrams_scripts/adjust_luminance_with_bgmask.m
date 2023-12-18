%% luminance control images
%  in case of questions please submit to the authour:
%  yeo-jin.yi.15@ucl.ac.uk

%% initialise

clear;clc
path = 'F:/ENGRAMS/emo/images/originals/scenes/';
path_save = 'F:/ENGRAMS/emo/images/scenes/';
cd(path)
dirData = dir('*.png');
fileNames = {dirData.name};
fileNames(ismember(fileNames,{'.','..','._.DS_Store','._*','.DS_Store'})) = [];
isHiddenFile = arrayfun(@(f) startsWith(f, '._'), fileNames);  % remove the bullshit from mac os
filteredFiles = fileNames(~isHiddenFile);  
fileNames = filteredFiles;

pictureOrigSize = [600 800]; % [height,width]
 
%% control

for i = 1:length(fileNames)
    dat = imread([path fileNames{i}]);
     
    R = im2double(squeeze(dat(:,:,1)))+.000001;% red, prevent values from getting zero
    G = im2double(squeeze(dat(:,:,2)))+.000001;% green
    B = im2double(squeeze(dat(:,:,3)))+.000001;% blue
    
    % make a background mask which is this shade of grey
    greyColor = [129, 129, 129];
    greyColorRep = repmat(reshape(greyColor, [1, 1, 3]), size(dat, 1), size(dat, 2));
    exclusionMask = all(dat == greyColorRep, 3); 

    rfact = 0.2126; gfact = 0.7152; bfact = 0.0722;
    Y(i,1) = rfact*geomean(R(:)) + gfact*geomean(G(:)) + bfact*geomean(B(:));
    cutval = .5;
    ratios = mean(log(R(:))) + mean(log(G(:))) + mean(log(B(:)));
     
    offset = .4;
     
    % adjust individual colors
    diff = mean(log(R(:)+offset))- log(cutval);
    newvals = exp(log(R(:)+offset)-diff);
    Rn = reshape(newvals,size(R));clear diff newvals
    diff = mean(log(G(:)+offset))- log(cutval);
    newvals = exp(log(G(:)+offset)-diff);
    Gn = reshape(newvals,size(G));clear diff newvals
    diff = mean(log(B(:)+offset))- log(cutval);
    newvals = exp(log(B(:)+offset)-diff);
    Bn = reshape(newvals,size(B));clear diff newvals

    % apply luminance adjustment except for the excluded grey color
    Rn(exclusionMask) = R(exclusionMask);
    Gn(exclusionMask) = G(exclusionMask);
    Bn(exclusionMask) = B(exclusionMask);
     
    Y(i,2) = rfact*geomean(Rn(:)) + gfact*geomean(Gn(:)) + bfact*geomean(Bn(:));
    dat2 = zeros(size(dat));
    dat2(:,:,1) = Rn;dat2(:,:,2) = Gn;dat2(:,:,3) = Bn;
    
    % display and save the image
    close all
    F = gcf;
    set(F,'Position', [100 300 768 576]);
    set(gca, 'position', [0 0 1 1], 'visible', 'off');
    imshow(dat2); pause(.5)
    G = getframe(gcf);
    
    % resize the image to whatever the size you want
    resizedImage = imresize(dat2, pictureOrigSize);

    % save the resized image
    imwrite(resizedImage, [path_save fileNames{i}]);
end

figure, plot(Y);
