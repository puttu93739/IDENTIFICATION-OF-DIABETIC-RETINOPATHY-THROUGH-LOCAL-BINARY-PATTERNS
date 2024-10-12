warning ('off')
clc;
clear all;
close all;

%%% PREPROCESSING%%%%
title('original image');
%E-OPTHA DATABASE%% (M X N X O)
%%% READING INPUT IMAGE%%%
[filename,pathname]=uigetfile( {'*.png'; '*.bmp';'*.tif';'*.jpg'});
I=imread([pathname filename]);
figure,imshow(I,[]);

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

