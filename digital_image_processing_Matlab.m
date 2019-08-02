function [ output_args ] = digital_image_processing_Matlab()
%Write your own code here
%When you have completed one task, you may comment output by 
%placing "%" in the beginning of row (if the result is not 
%needed later)
%if you get a lot of prints in Matlab console (slows your code), put ; to the end of line

%readfile
close all ;
clc;

imfinfo('lena.tif');
RGBimage=imread('lena.tif');
imshow(RGBimage)

%show the image and its information
figure(1),imshow(RGBimage),impixelinfo %xy coordinate and RGB value
impixel(RGBimage,100, 100) %RGB value on specific xycoordinate

%turn gray
figure(2),imshow(RGBimage),impixelinfo
figure(3),imshow(RGBimage(:,:,1)) %show the gray level image
grayimage=rgb2gray(RGBimage); %turn the image into gray level
figure(4),imshow(grayimage)

%gray level
grayimage16=grayslice(grayimage, 16);% gray scale 16
imshow(grayimage16)
grayimage4=grayslice(grayimage, 4); % gray scale 4
imshow(grayimage4)
figure(5), imhist(grayimage16),axis([-50 300 0 max(imhist(grayimage16))])
figure(6), imhist(grayimage4),axis([-50 300 0 max(imhist(grayimage4))])% show the distribution of image
figure(7),imshow(histeq(grayimage16))% make the gray level image evenly distributed
figure(8),imshow(histeq(grayimage4))


mycolormap=[0.5 0 1; 1 0.5 1;0 1 0.5; 1 0.5 0.5];
figure(9),imshow(grayimage4,mycolormap),impixelinfo %turn gray scale to another color

%binary image
binimage=grayimage<150;
figure(10),imshow(binimage),impixelinfo

pause

%noise image
noiseimage=imnoise(grayimage,'salt & pepper');% give salt and peper to the image
imshow(noiseimage)

d_grayimage=im2double(grayimage);% turn the image into double
imshow(d_grayimage)
newimage=1-d_grayimage;
imshow(newimage)
imshow(imsubtract(grayimage,128)) %minus the image
imshow(imadd(grayimage,128))% add
imshow(imdivide(grayimage,2))%divide
imshow(immultiply(RGBimage,2))%multiply

pause

%mask image
maskimage=imread('lena_mask.tif');
imshow(maskimage)
binimage2=maskimage<150;
figure(3),imshow(binimage2),impixelinfo;
imshow(immultiply(grayimage,binimage2))% show the image in the mask

pause

%convolution
a=ones(3,3)/9
figure(4),imshow(conv2(double(grayimage),double(a))/100),impixelinfo;
b=ones(5,5)/25;
figure(5),imshow(conv2(double(grayimage),double(b))/100),impixelinfo;
c=ones(9,9)/81;
figure(6),imshow(conv2(double(grayimage),double(c))/100),impixelinfo;

pause

%Laplacian
mask=fspecial('laplacian');
mask = [0.1667 0.6667 0.1667; 0.6667 -3.3333 0.6667; 0.1667 0.6667 0.1667];
figure(7),imshow(conv2(double(grayimage),double(mask))),impixelinfo;
figure(8),imshow(conv2(double(grayimage),double(mask))/100),impixelinfo;

p=fspecial('average');% give a gaussian matrix
average_image=conv2(double(grayimage),double(p));
average_image= average_image(2:513,2:513,:);
figure(9),imshow(average_image),impixelinfo;

figure(10),imshow(double(grayimage)-average_image/1.5),impixelinfo;
figure(11),imshow((double(grayimage)-average_image/1.5)/70),impixelinfo;

pause

%fourier
f_image=fftshift(fft2(grayimage));
figure(12),imshow(mat2gray(log(1+abs(f_image)))),impixelinfo %fast fourier transformation

[x,y]=meshgrid(-256:255,-256:255);
z=sqrt(x.^2+y.^2);
c=z<15;
lowpass_f_image=f_image.*c;
figure(13),imshow(mat2gray(log(1+abs(lowpass_f_image)))),impixelinfo % low pass filter on fourier transformation
lowpass_if_image=ifft2(ifftshift(lowpass_f_image));%turn back
figure(14),imshow(uint8(lowpass_if_image)),impixelinfo

c=z<50;
lowpass_f_image=f_image.*c;
figure(15),imshow(mat2gray(log(1+abs(lowpass_f_image)))),impixelinfo
lowpass_if_image=ifft2(ifftshift(lowpass_f_image));
figure(16),imshow(uint8(lowpass_if_image)),impixelinfo

gaussianfilter=mat2gray(fspecial('gaussian',512,20)); %gaussian filter on fourier transformation
lowpass_gf_image=f_image.*gaussianfilter;
figure(17),imshow(lowpass_gf_image),impixelinfo
lowpass_igf_image=ifft2(ifftshift(lowpass_gf_image));
figure(18),imshow(lowpass_igf_image/100),impixelinfo

pause

c=z>15;
highpass_f_image=f_image.*c; %high pass filter on fourier transformation
figure(19),imshow(mat2gray(log(1+abs(highpass_f_image)))),impixelinfo
highpass_if_image=ifft2(ifftshift(highpass_f_image));
figure(20),imshow(uint8(highpass_if_image)),impixelinfo

c=z>50;
highpass_f_image=f_image.*c;
figure(21),imshow(mat2gray(log(1+abs(highpass_f_image)))),impixelinfo
highpass_if_image=ifft2(ifftshift(highpass_f_image));
figure(22),imshow(uint8(highpass_if_image)),impixelinfo

pause

gas_mask=1-mask;
figure(23),imshow(conv2(double(grayimage),double(gas_mask))),impixelinfo
figure(24),imshow(conv2(double(grayimage),double(gas_mask))/500),impixelinfo

s=size(grayimage);
[x,y]=meshgrid(1:s(1),1:s(2)); %give noise to the image and fourier
p=sin(x/7+y/6)+1;
periodic_noiseimage=(im2double(grayimage)+p/2)/2;
fp=fftshift(fft2(periodic_noiseimage));
figure(25),imshow(mat2gray(log(1+abs(fp)))),impixelinfo

f_periodic_noiseimage=fp;
z=sqrt((x-256).^2+(y-256).^2);
bandreject_mask=(z < 0.45|z>0.5); %band reject and fourier
band_rejected_image= f_periodic_noiseimage.*bandreject_mask;
figure(26),imshow(band_rejected_image),impixelinfo
figure(27),imshow(ifft2(ifftshift(band_rejected_image))),impixelinfo

band_stop_filter=ones(512,512);% band stop filter on fourier
band_stop_filter(244:246,:)=0;
band_stop_filter(269:271,:)=0;
band_stop_filter(:,244:246)=0;
band_stop_filter(:,269:271)=0;
band_stopped_f_periodic_noiseimage= f_periodic_noiseimage.* band_stop_filter;
figure(28),imshow(band_stopped_f_periodic_noiseimage),impixelinfo
figure(29),imshow(ifft2(ifftshift(band_stopped_f_periodic_noiseimage))),impixelinfo


pause


close all
end

