function imageOut(max_area,area_Yeast,area_Bacteria,count, flag,ImageYeast, Ya_layer, ImageBacteria, Bacteria_layer, ImageBacteria_temp2,ImageYeast_temp2, ImageYeast_Fluid)
n = 512;
Ya_layer = imresize(Ya_layer,[n, nan],'nearest');
Bacteria_layer = imresize(Bacteria_layer,[n, nan],'nearest');
ImageBacteria_temp2 = imresize(ImageBacteria_temp2,[n, nan],'nearest');
ImageYeast_temp2 = imresize(ImageYeast_temp2,[n, nan],'nearest');
ImageYeast_Fluid = imresize(ImageYeast_Fluid,[n, nan],'nearest');
ImageYeast = imresize(ImageYeast,[n, nan],'nearest');
ImageBacteria = imresize(ImageBacteria,[n, nan],'nearest');
colorImg = ones([n,n,3]);

loc = strcat(pwd,'\Frames');
%% Output images

hYA = 0.13;
sYA = 0.349;
map_YA = hsv2rgb([hYA,sYA,0.5;    hYA,sYA,0.55;    hYA,sYA,0.6;    hYA,sYA,0.65;    hYA,sYA,0.7;
    hYA,sYA,0.75;    hYA,sYA,0.8;    hYA,sYA,0.85;    hYA,sYA,0.9;
    hYA,sYA,0.95;    hYA,sYA,1;    ]);

hPA = 0.3;
sPA = 0.2706;
map_PA = hsv2rgb([hPA,sPA,0;    hPA,sPA,0.1;    hPA,sPA,0.2;
    hPA,sPA,0.3;    hPA,sPA,0.4;    hPA,sPA,0.5;    hPA,sPA,0.6;
    hPA,sPA,0.7;    hPA,sPA,0.8;    hPA,sPA,0.9;    hPA,sPA,1;    ]);

%YA space hsv(35,89,100)
Rind_blk = find(Ya_layer == 1);
colorImg(Rind_blk) = hYA;
colorImg(Rind_blk+(n*n)) = sYA;
%     colorImg(Rind_blk+20000) = 2.*(ImageYeast_temp2(Rind_blk)./max(max(ImageYeast_temp2)));
colorImg(Rind_blk+(2*n*n)) = 1 - (ImageYeast(Rind_blk)./(max_area/area_Yeast)).*0.5;
clear Rind_blk

Rind_blk = find(Ya_layer == 0);
colorImg(Rind_blk) = 1;% hYA-0.03;
colorImg(Rind_blk+(n*n)) = 0;%sYA+0.5; %colorImg(Rind_blk+10000) = 0.1;
colorImg(Rind_blk+(2*n*n)) = 1;%0.9; %colorImg(Rind_blk+20000) = 1
colorImg2 = colorImg;
imgg_YA = hsv2rgb(colorImg2);
clear Rind_blk colorImg2

%blank space hsv(0,0,100)
Rind_blk = find(Ya_layer == 0 & Bacteria_layer == 0);
colorImg(Rind_blk) = 1;% hYA-0.03;
colorImg(Rind_blk+(n*n)) = 0;%sYA+0.5;
colorImg(Rind_blk+(2*n*n)) = 1;%0.9;
clear Rind_blk

%PA space hsvhsv(159,39,38)
Rind_blk = find(Bacteria_layer == 1 );
colorImg(Rind_blk) = hPA;%0.8
colorImg(Rind_blk+(n*n)) = sPA;%0.8
% colorImg(Rind_blk+20000) = (colorImg(Rind_blk+20000)+(ImageSwarm_temp2(Rind_blk)./max(max(ImageSwarm_temp2))));
Temp_pa_imag = 1 - (ImageBacteria(Rind_blk)/(max_area/area_Bacteria)).*0.5;
Temp_ya_imag = 0 ;%+ ((ImageYeast_temp2(Rind_blk)./max(max(ImageYeast_temp2)))./2);
colorImg(Rind_blk+(2*n*n)) = Temp_ya_imag + Temp_pa_imag;

temp_imageReszed(:,:,1) = imresize(colorImg(:,:,1),[n, nan],'nearest');
temp_imageReszed(:,:,2) = imresize(colorImg(:,:,2),[n, nan],'nearest');
temp_imageReszed(:,:,3) = imresize(colorImg(:,:,3),[n, nan],'nearest');
imgg = hsv2rgb(temp_imageReszed);
map = colormap("sky");
img = ImageYeast_Fluid;
bg = [0 0 0];
valid = isfinite(img);
bounds = [min(img(:)) max(img(:))];
% bounds = [0 4e-3];%[min(img(:)) max(img(:))];%[0 4e-3];%[3.7e-3,4.3e-3];
colorbar;
n = size(map, 1);
inds = (img(valid)-bounds(1))/(bounds(end)-bounds(1))*(n-1);
inds = floor(min(max(inds, 0), n-1))+1;
dim = size(img);
r = ones(dim)*bg(1); r(valid) = map(inds, 1);
g = ones(dim)*bg(2); g(valid) = map(inds, 2);
b = ones(dim)*bg(3); b(valid) = map(inds, 3);
dim2 = [dim(1:2) 3 dim(3:end)];
rgb = zeros(dim2, class(map));
rgb(:,:,1,:) = r;
rgb(:,:,2,:) = g;
rgb(:,:,3,:) = b;
% imshow(rgb,[])
% imshow(rgb)

loc1 = strcat(loc,'\YeastLayer\',num2str(flag),'\',num2str(count),'.jpg');
loc2 = strcat(loc,'\CoExit\',num2str(flag),'\',num2str(count),'.jpg');
loc3 = strcat(loc,'\FluidLayer\',num2str(flag),'\',num2str(count),'.jpg');
%% export
imwrite(imgg_YA,loc1);
imwrite(imgg,loc2);
imwrite(rgb,loc3)
end
