warning ('off')
clc;
clear all;
close all;

%%% PREPROCESSING%%%%

%E-OPTHA DATABASE%% (M X N X O)
%%% READING INPUT IMAGE%%%
[filename,pathname]=uigetfile( {'*.png'; '*.bmp';'*.tif';'*.jpg'});
I=imread([pathname filename]);
figure,imshow(I,[]);
title('original image');

%% RESCALING%%%%
%%%%%%%% take bicubic interpolation%%%%%%%%
%%%%%%%%resize the image size%%%%%%%%%%%%%
RI=imresize(I,[512,512],'bicubic');
figure,imshow(RI,[]);
title('resized input image');
impixelinfo;

%% CHANNEL SEPARATION %%%
red=RI(:,:,1);
green=RI(:,:,2);
blue=RI(:,:,3);

figure,imshow(red,[]);
title('red channel image');
impixelinfo;

figure,imshow(green,[]);
title('green channel image');
impixelinfo;

figure,imshow(blue,[]);
title('blue channel image');
impixelinfo;


%% MEDIAN FILTER %%%
%%% median filter-red channel%%

filtred=medfilt2(red,[3 3]);
[rr,rc]=size(filtred);

for ri=1:rr
    for rj=1:rc
        if(filtred(ri,rj)<4)
            filtred(ri,rj)=0;
        else
        end
    end
end


figure,imshow(filtred,[]);
title('red channel median filtered image');
impixelinfo;
[rr rc rd]=size(filtred);

%% median filter-green channel%%
filtgreen=medfilt2(green,[3 3]);

[gr gc]=size(filtgreen);

for gi=1:gr
    for gj=1:gc
        if(filtgreen(gi,gj)<4)
            filtgreen(gi,gj)=0;
        else
        end
    end
end


figure,imshow(filtgreen,[]);
title('green channel median filtered image');
impixelinfo;

%% median filter-blue channel%%
filtblue=medfilt2(blue,[3 3]);

[br bc]=size(filtblue);

for bi=1:br
    for bj=1:bc
        if(filtblue(bi,bj)<4)
            filtblue(bi,bj)=0;
        else
        end
    end
end

figure,imshow(filtblue,[]);
title('blue channel median filtered image');
impixelinfo;


%% %MASKING%%%

%%%%OPTIC DISK AND VESSEL REMOVAL%%%%
%%% OPTIC DISK AND BACKGROUND REMOVAL FOR RED CHANNEL%%%

%%
% %%%% vessel mask%%%%
%Extract Blood Vessels
Threshold = 5;
rbloodVesselsr = VesselExtract(filtred, Threshold);
figure,imshow(rbloodVesselsr);
title('Extracted Blood Vessels-RED');

bloodVesselsr=im2bw(rbloodVesselsr);
figure,imshow(bloodVesselsr);
title('Extracted Blood Vessels-RED BW');

se=strel('disk',1);
bloodVesselsr=imdilate(bloodVesselsr,se);
figure,imshow(bloodVesselsr);
title('holesfilled image');

bloodVesselsri=imcomplement(bloodVesselsr);
figure,imshow(bloodVesselsri);
ROIvr(~bloodVesselsri)=0;
title('Vessel mask-RED');

%%% masking with the input image
ROIvr=filtred;
figure,imshow(ROIvr);
title('vessel removed image-red channel');
impixelinfo;

%% DISC REMOVAL PROCESS FOR RED CHANNEL%%%

rmu = mean2(red); % Computing the mean
%%%%%%%%%%% Calculating the mean subtraction from the input image%%%
rmeansub=red-rmu; % Computing the difference image for each image Ai = Ti - m
figure,imshow(rmeansub);
title('mean subtracted image-red channel');
impixelinfo;

rROI1=rmeansub(:,:,1)<100; 
figure,imshow(rROI1);
title('ROI1 MASK');

%%% masking with the input image
% ROIred=filtred;
ROIred=ROIvr;
ROIred(~rROI1)=0;
figure,imshow(ROIred);
title('optic disk removed image-red channel');
impixelinfo;

%% % OPTIC DISK AND BLOOD VESSEL REMOVAL FOR GREEN CHANNEL%%%

%Extract Blood Vessels
Threshold = 6;
gbloodVesselsg = VesselExtract(filtgreen, Threshold);
figure,imshow(gbloodVesselsg);
title('Extracted Blood Vessels-GREEN');

bloodVesselsg=im2bw(gbloodVesselsg);
figure,imshow(bloodVesselsg);
title('Extracted Blood Vessels-GREEN BW');

se=strel('disk',1);
bloodVesselsg=imdilate(bloodVesselsg,se);
figure,imshow(bloodVesselsg);
title('holes filled image');

%%% masking with the input image
bloodVesselsgi=imcomplement(bloodVesselsg);
figure,imshow(bloodVesselsgi);
title('Vessel mask-GREEN');

ROIvg=filtgreen;
ROIvg(~bloodVesselsgi)=0;
figure,imshow(ROIvg);
title('vessel removed image-green channel');
impixelinfo;

%% DISC REMOVAL PROCESS FOR Green CHANNEL%%%
gmu = mean2(green); % Computing the mean
%%%%%%%%%%% Calculating the mean subtraction from the input image%%%
gmeansub=green-gmu; % Computing the difference image for each image Ai = Ti - m
figure,imshow(gmeansub);
title('mean subtracted image--green channel');
impixelinfo;

ROI2=gmeansub(:,:,1)<60; 
figure,imshow(ROI2);
title('ROI1 MASK-GREEN');

ROIgreen=ROIvg;
ROIgreen(~ROI2)=0;
figure,imshow(ROIgreen);
title('optic disk removed image-green channel');
impixelinfo;

%% % OPTIC DISK AND VESSEL BACKGROUND REMOVAL FOR BLUE CHANNEL%%%
%Extract Blood Vessels
Threshold=10;
bbloodVesselsb = VesselExtract(filtblue, Threshold);
figure,imshow(bbloodVesselsb);
title('Extracted Blood Vessels');

bloodVesselsb=im2bw(bbloodVesselsb);
figure,imshow(bloodVesselsb);
title('Extracted Blood Vessels');

se=strel('disk',1);
bloodVesselsb=imdilate(bloodVesselsb,se);
figure,imshow(bloodVesselsb);
title('holes filled image');

%%% masking with the input image
bloodVesselsbi=imcomplement(bloodVesselsb);
figure,imshow(bloodVesselsbi);
title('Vessel mask');

ROIvb=filtblue;
ROIvb(~bloodVesselsbi)=0;
figure,imshow(ROIvb);
title('vessel removed image-blue channel');
impixelinfo;

%% DISC REMOVAL PROCESS FOR Blue CHANNEL%%%
%%% Background subtraction%%%
bmu = mean2(blue); % Computing the mean
%%%%%%%%%%% Calculating the mean subtraction from the input image%%%
bmeansub=blue-bmu; % Computing the difference image for each image Ai = Ti - m
figure,imshow(bmeansub);
title('mean subtracted image--blue channel');
impixelinfo;

bROI1=bmeansub(:,:,1)<30; 
figure,imshow(bROI1);
title('ROI1 MASK');

ROIblue=ROIvb;
ROIblue(~bROI1)=0;
figure,imshow(ROIblue);
title('optic disk removed image-blue channel');
impixelinfo;

%% %% FEATURE EXTRACTION %%%%%
%%% TEXTURE DESCRIPTORS%%%%
%%%%LOCAL BINARY PATTERNS%%%%
%%%%LBP WITH DIFFERENT RADIUS%%%%
%%% NEIGHBOR P=8,RADIUS=1 %%%%
%%% RED CHANNEL--RADIUS 1%%%
filtred=ROIred;
nFiltSize=8;
nFiltRadius=1;
filtR1=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPR1=efficientLBP(filtred, 'filtR1', filtR1, 'isRotInv', false, 'isChanWiseRot', false);
effLBPR1=imcomplement(effLBPR1);
figure,imshow(effLBPR1);
title('LBP feature extracted Red channel image-1');
impixelinfo;

%% GREEN CHANNEL--RADIUS 1%%%
filtgreen=ROIgreen;
nFiltSize=8;
nFiltRadius=1;
filtG1=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPG1=efficientLBP(filtgreen, 'filtG1', filtG1, 'isRotInv', false, 'isChanWiseRot', false);
effLBPG1=imcomplement(effLBPG1);
figure,imshow(effLBPG1);
title('LBP feature extracted Green channel image-1');

%% BLUE CHANNEL--RADIUS 1%%%
filtblue=ROIblue;
nFiltSize=8;
nFiltRadius=1;
filtB1=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPB1=efficientLBP(filtblue, 'filtB1', filtB1, 'isRotInv', false, 'isChanWiseRot', false);
effLBPB1=imcomplement(effLBPB1);
figure,imshow(effLBPB1);
title('LBP feature extracted Blue channel image-1');

%% NEIGHBOR P=8,RADIUS=2 %%%%
%%% RED CHANNEL--RADIUS 2%%%
filtred=ROIred;
nFiltSize=8;
nFiltRadius=2;
filtR2=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPR2=efficientLBP(filtred, 'filtR2', filtR2, 'isRotInv', false, 'isChanWiseRot', false);
effLBPR2=imcomplement(effLBPR2);
figure,imshow(effLBPR2);
title('LBP feature extracted Red channel image-2');

%% GREEN CHANNEL--RADIUS 2%%%
filtgreen=ROIgreen;
nFiltSize=8;
nFiltRadius=2;
filtG2=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPG2=efficientLBP(filtgreen, 'filtG2', filtG2, 'isRotInv', false, 'isChanWiseRot', false);
effLBPG2=imcomplement(effLBPG2);
figure,imshow(effLBPG2);
title('LBP feature extracted Green channel image-2');

%% BLUE CHANNEL--RADIUS 2%%%
filtblue=ROIblue;
nFiltSize=8;
nFiltRadius=2;
filtB2=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPB2=efficientLBP(filtblue, 'filtB2', filtB2, 'isRotInv', false, 'isChanWiseRot', false);
effLBPB2=imcomplement(effLBPB2);
figure,imshow(effLBPB2);
title('LBP feature extracted Blue channel image-2');

%% % NEIGHBOR P=8,RADIUS=3 %%%%
%%% RED CHANNEL--RADIUS 3 %%%
filtred=ROIred;
nFiltSize=8;
nFiltRadius=3;
filtR3=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPR3=efficientLBP(filtred, 'filtR3', filtR3, 'isRotInv', false, 'isChanWiseRot', false);
effLBPR3=imcomplement(effLBPR3);
figure,imshow(effLBPR3);
title('LBP feature extracted Red channel image-3');

%% GREEN CHANNEL--RADIUS 3%%%
filtgreen=ROIgreen;
nFiltSize=8;
nFiltRadius=3;
filtG3=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPG3=efficientLBP(filtgreen, 'filtG3', filtG3, 'isRotInv', false, 'isChanWiseRot', false);
effLBPG3=imcomplement(effLBPG3);
figure,imshow(effLBPG3);
title('LBP feature extracted Green channel image-3');

%% BLUE CHANNEL--RADIUS 3%%%
filtblue=ROIblue;
nFiltSize=8;
nFiltRadius=3;
filtB3=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPB3=efficientLBP(filtblue, 'filtB3', filtB3, 'isRotInv', false, 'isChanWiseRot', false);
effLBPB3=imcomplement(effLBPB3);
figure,imshow(effLBPB3);
title('LBP feature extracted Blue channel image-3');

%% % NEIGHBOR P=8,RADIUS=5 %%%%
%%% RED CHANNEL--RADIUS 5 %%%
filtred=ROIred;
nFiltSize=8;
nFiltRadius=5;
filtR5=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPR5=efficientLBP(filtred, 'filtR5', filtR5, 'isRotInv', false, 'isChanWiseRot', false);
effLBPR5=imcomplement(effLBPR5);
figure,imshow(effLBPR5);
title('LBP feature extracted Red channel image-5');

%% GREEN CHANNEL--RADIUS 5%%%
filtgreen=ROIgreen;
nFiltSize=8;
nFiltRadius=5;
filtG5=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPG5=efficientLBP(filtgreen, 'filtG5', filtG5, 'isRotInv', false, 'isChanWiseRot', false);
effLBPG5=imcomplement(effLBPG5);
figure,imshow(effLBPG5);
title('LBP feature extracted Green channel image-5');

%% BLUE CHANNEL--RADIUS 5%%%
filtblue=ROIblue;
nFiltSize=8;
nFiltRadius=5;
filtB5=generateRadialFilterLBP(nFiltSize, nFiltRadius);
effLBPB5=efficientLBP(filtblue, 'filtB5', filtB5, 'isRotInv', false, 'isChanWiseRot', false);
effLBPB5=imcomplement(effLBPB5);
figure,imshow(effLBPB5);
title('LBP feature extracted Blue channel image-5');

%% LBP HISTOGRAM%%%

%% RADIUS-1 NEIGHBOR-8%%%
Rlbpout1=lbp(filtred,1,8);
Glbpout1=lbp(filtgreen,1,8);
Blbpout1=lbp(filtblue,1,8);

%% RADIUS-2 NEIGHBOR-8%%%
Rlbpout2=lbp(filtred,2,8);
Glbpout2=lbp(filtgreen,2,8);
Blbpout2=lbp(filtblue,2,8);

%% RADIUS-3 NEIGHBOR-8%%%
Rlbpout3=lbp(filtred,3,8);
Glbpout3=lbp(filtgreen,3,8);
Blbpout3=lbp(filtblue,3,8);

%% RADIUS-5 NEIGHBOR-8%%%
Rlbpout5=lbp(filtred,5,8);
Glbpout5=lbp(filtgreen,5,8);
Blbpout5=lbp(filtblue,5,8);

%% %% VARIANCE (VAR) %%%%
%%%  NEIGHBOR P=8,RADIUS=1 %%%%
varoutR1=varn(filtred,1,8);
varoutR1=imcomplement(varoutR1);
figure,imshow(varoutR1, []);
title('Variance Feature Image-Red Channel  Radius-1');
hvaroutR1=imhist(varoutR1);

varoutG1=varn(filtgreen,1,8);
varoutG1=imcomplement(varoutG1);
figure,imshow(varoutG1,[]);
title('Variance Feature Image-Green Channel  Radius-1');
hvaroutG1=imhist(varoutG1);

varoutB1=varn(filtblue,1,8);
varoutB1=imcomplement(varoutB1);
figure,imshow(varoutB1,[]);
title('Variance Feature Image-Blue Channel  Radius-1');
hvaroutB1=imhist(varoutB1);

%% %% VARIANCE (VAR) %%%%
%%%  NEIGHBOR P=8,RADIUS=2 %%%%
varoutR2=varn(filtred,2,8);
varoutR2=imcomplement(varoutR2);
figure,imshow(varoutR1, []);
title('Variance Feature Image-Red Channel  Radius-2');
hvaroutR2=imhist(varoutR2);

varoutG2=varn(filtgreen,2,8);
varoutG2=imcomplement(varoutG2);
figure,imshow(varoutG1,[]);
title('Variance Feature Image-Green Channel  Radius-2');
hvaroutG2=imhist(varoutG2);

varoutB2=varn(filtblue,2,8);
varoutB2=imcomplement(varoutB2);
figure,imshow(varoutB2,[]);
title('Variance Feature Image-Blue Channel  Radius-2');
hvaroutB2=imhist(varoutB2);

%% %% VARIANCE (VAR) %%%%
%%%  NEIGHBOR P=8,RADIUS=3 %%%%
varoutR3=varn(filtred,3,8);
varoutR3=imcomplement(varoutR3);
figure,imshow(varoutR3, []);
title('Variance Feature Image-Red Channel Radius-3');
hvaroutR3=imhist(varoutR3);

varoutG3=varn(filtgreen,3,8);
varoutG3=imcomplement(varoutG3);
figure,imshow(varoutG3,[]);
title('Variance Feature Image-Green Channel  Radius-3');
hvaroutG3=imhist(varoutG3);

varoutB3=varn(filtblue,3,8);
varoutB3=imcomplement(varoutB3);
figure,imshow(varoutB3,[]);
title('Variance Feature Image-Blue Channel  Radius-3');
hvaroutB3=imhist(varoutB3);

%% %% VARIANCE (VAR) %%%%
%%%  NEIGHBOR P=8,RADIUS=5 %%%%
varoutR5=varn(filtred,5,8);
varoutR5=imcomplement(varoutR5);
figure,imshow(varoutR5, []);
title('Variance Feature Image-Red Channel Radius-5');
hvaroutR5=imhist(varoutR5);

varoutG5=varn(filtgreen,5,8);
varoutG5=imcomplement(varoutG5);
figure,imshow(varoutG5,[]);
title('Variance Feature Image-Green Channel  Radius-5');
hvaroutG5=imhist(varoutG5);

varoutB5=varn(filtblue,5,8);
varoutB5=imcomplement(varoutB5);
figure,imshow(varoutB5,[]);
title('Variance Feature Image-Blue Channel  Radius-5');
hvaroutB5=imhist(varoutB5);

%% FEATURE measurement%%%
[mu1,sd1,med1,ent1,skew1,kurt1]=texturefeature(Rlbpout1);
[mu2,sd2,med2,ent2,skew2,kurt2]=texturefeature(Glbpout1);
[mu3,sd3,med3,ent3,skew3,kurt3]=texturefeature(Blbpout1);

[mu4,sd4,med4,ent4,skew4,kurt4]=texturefeature(Rlbpout2);
[mu5,sd5,med5,ent5,skew5,kurt5]=texturefeature(Glbpout2);
[mu6,sd6,med6,ent6,skew6,kurt6]=texturefeature(Blbpout2);

[mu7,sd7,med7,ent7,skew7,kurt7]=texturefeature(Rlbpout3);
[mu8,sd8,med8,ent8,skew8,kurt8]=texturefeature(Glbpout3);
[mu9,sd9,med9,ent9,skew9,kurt9]=texturefeature(Blbpout3);

[mu10,sd10,med10,ent10,skew10,kurt10]=texturefeature(Rlbpout5);
[mu11,sd11,med11,ent11,skew11,kurt11]=texturefeature(Glbpout5);
[mu12,sd12,med12,ent12,skew12,kurt12]=texturefeature(Blbpout5);

[mu13,sd13,med13,ent13,skew13,kurt13]=texturefeature(hvaroutR1);
[mu14,sd14,med14,ent14,skew14,kurt14]=texturefeature(hvaroutG1);
[mu15,sd15,med15,ent15,skew15,kurt15]=texturefeature(hvaroutB1);

[mu16,sd16,med16,ent16,skew16,kurt16]=texturefeature(hvaroutR2);
[mu17,sd17,med17,ent17,skew17,kurt17]=texturefeature(hvaroutG2);
[mu18,sd18,med18,ent18,skew18,kurt18]=texturefeature(hvaroutB2);

[mu19,sd19,med19,ent19,skew19,kurt19]=texturefeature(hvaroutR3);
[mu20,sd20,med20,ent20,skew20,kurt20]=texturefeature(hvaroutG3);
[mu21,sd21,med21,ent21,skew21,kurt21]=texturefeature(hvaroutB3);

[mu22,sd22,med22,ent22,skew22,kurt22]=texturefeature(hvaroutR5);
[mu23,sd23,med23,ent23,skew23,kurt23]=texturefeature(hvaroutG5);
[mu24,sd24,med24,ent24,skew24,kurt24]=texturefeature(hvaroutB5);

%% %%% CLASSIFICATION%%%%
%%% Classification Using SVM

totalImages1=16;    
folder_name1='G:\final project reports\CODING-100%-full\images\TRAINING\e_ophtha_MA\healthy';
jpgImagesDir1 = fullfile(folder_name1, '*.jpg');
num_of_jpg_images1 = numel( dir(jpgImagesDir1));
jpg_files1 = dir(jpgImagesDir1);
% if ( ~isempty( jpg_files )
jpg_counter1 = 0;

 for k = 1:totalImages1
        
        if ( (num_of_jpg_images1 - jpg_counter1) > 0)
            imgInfoJPG1 = imfinfo( fullfile(folder_name1,jpg_files1(jpg_counter1+1).name ) );
            if ( strcmp( lower(imgInfoJPG1.Format), 'jpg') == 1 )
                % read images
                sprintf('%s \n', jpg_files1(jpg_counter1+1).name)
                % extract features
                image1 = imread( fullfile(folder_name1, jpg_files1(jpg_counter1+1).name ) );
                [pathstr1, name1, ext1] = fileparts( fullfile(folder_name1, jpg_files1(jpg_counter1+1).name ) );
                 image1 = imresize(image1, [512 512]);
            end
            
            jpg_counter1 = jpg_counter1 + 1;
        end
 end       

totalImages2=8;
folder_name2='G:\final project reports\CODING-100%-full\images\TRAINING\e_ophtha_MA\MA';
jpgImagesDir2 = fullfile(folder_name2, '*.jpg');
num_of_jpg_images2 = numel( dir(jpgImagesDir2));
jpg_files2 = dir(jpgImagesDir2);
% if ( ~isempty( jpg_files )
jpg_counter2 = 0;

 for k = 1:totalImages2
        
        if ( (num_of_jpg_images2 - jpg_counter2) > 0)
            imgInfoJPG2 = imfinfo( fullfile(folder_name2,jpg_files2(jpg_counter2+1).name ) );
            if ( strcmp( lower(imgInfoJPG2.Format), 'jpg') == 1 )
                % read images
                sprintf('%s \n', jpg_files2(jpg_counter2+1).name)
                % extract features
                image2 = imread( fullfile(folder_name2, jpg_files2(jpg_counter2+1).name ) );
                [pathstr2, name2, ext2] = fileparts( fullfile(folder_name2, jpg_files2(jpg_counter2+1).name ) );
                 image2 = imresize(image2, [512 512]);
            end
            
            jpg_counter2 = jpg_counter2 + 1;
        end
 end     
 
%%%%  Load the training set
%%% Taking six features for test
%%% Mean Standard Deviation, Median, Entropy, skewness, Kurtosis%%
testfeat_eyedisease=[mu1 sd1 med1 ent1 skew1 kurt1];
Train_Feat=[1016,9600.6,5,0.8822,15.5750,246.8904; 1016,9338.5,4,0.9485,15.5055,245.3771]; 
Train_Label=[0;1];
% % Train the classifier 
svmStructDisease = svmtrain(Train_Feat,Train_Label);
%%% classify the retinal image
type_disease = svmclassify(svmStructDisease,testfeat_eyedisease);
%  Observe the results on the command window
disp(testfeat_eyedisease);
if type_disease == 0
    helpdlg('healthy image');
    disp(' healthy image');
elseif type_disease == 1
    helpdlg('MA affected image');
    disp('MA affected image');
end

