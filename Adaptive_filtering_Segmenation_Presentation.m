%% adding need java library
    javaaddpath 'ij.jar';
    javaaddpath 'mij.jar';
%% Clearing space
clear all;close all;clc;
%% Read the TIF file
[FileName,PathName] = uigetfile({'*.tif';'*.tiff';'*.png';'*.jpg';'*.bmp'},'Select the input tif file');
dname = PathName;
name{1} = FileName;
kid=[1]; % set as 1 or 2
fname=name{kid};
mkdir([dname filesep strrep(fname,'.tif','') '_Results'])
nametex=[dname filesep strrep(fname,'.tif','') '_Results' filesep strrep(fname,'.tif','')];

info = imfinfo([PathName fname]);
num_images = numel(info);
I = imread([PathName fname]);
I=I(:,:,1);


    
NCH=1; %'Is the staining of the border? 1 = Yes, 0 = No';
NR=1; %'WANT TO DO NORMALIZATION? 1 = Yes, 0 = No';
NCHmax=2; %'Minimum valid width of borders in pixels?';
MAXS=50000; % Maximum size of object in pixels
LTH=0.08;

SM=0; %% need smoothing? 1 = Yes, 0 = No;
if SM==1
h = fspecial('gaussian',[21 21],3);
I = imfilter(I,h,'replicate');
end

imwrite(uint16(65535*mat2gray(I)),[nametex 'Ori.png']);  

if NCH==1
    
     if NR==1
     I=uint16(65535*mat2gray(localnormalize(double(I),75,75)));
     end
     
     Img1=I;
        mean_s=NCHmax/2;  
        mins=0.51;
        maxs=max([mins NCHmax/2]);
        kernel_s=round(maxs*5);
    
  IM4=Img1; 
  fname1=strrep(fname,'.tif','');  
  Img=Img1;

        h = fspecial('disk', mean_s);
        Img2 = imfilter(Img,h,'replicate');

        Img3=imcomplement(Img2);
        Img8=imcomplement(Img);
        max3=max(Img3(:));
        min3=min(Img3(:));

        type=16;
        types=linspace(mins,maxs,type);
              Chd= uint16(round(type-1)*mat2gray(Img3))+1;
              Imgn=[];
                   for k=1:type
                  
                      h = fspecial('gaussian',[kernel_s kernel_s],types(k));
                                       Imgn(:,:,k) = imfilter(Img8,h,'replicate');
                   end
              adapt_image=zeros(size(Img8));
              
              for k=1:type
                        temp=Imgn(:,:,k);
                      adapt_image(Chd==k)=temp(Chd==k);
                                  
              end

        Img1=uint16(65535*mat2gray(adapt_image));
        imwrite(Img1,[nametex 'Adaptive_gaussian_filtered' '.png']);
                 
       Im2=Img1;
       Im4=uint16(65535*mat2gray(IM4));
       idx=1; 
       tol(idx)=std(double(Im2(:)))/sqrt(3);
       
    MIJ.createImage(Im2);
    MIJ.run('Find Maxima...', sprintf('noise=%d output=[Segmented Particles]',tol(idx)));
    I = MIJ.getCurrentImage;
    
    MIJ.closeAllWindows;
    Img11=uint16(65535*mat2gray(I));
    ImgOut11=uint16(plot_maxima2(Im4,imcomplement(Img11)));

                imwrite(uint16(65535*mat2gray(ImgOut11)),[nametex 'Segmenation.png']);
                imwrite(uint16(65535*mat2gray(Img11)),[nametex 'Border.png']);
else
            
IO=I;        
gradmagk=[];
hG = fspecial('gaussian',[15 15],0.1);
for kin=1:1    
    I(:,:,kin)=imfilter(I(:,:,kin),hG,'replicate');
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I(:,:,kin)), hy, 'replicate');
Ix = imfilter(double(I(:,:,kin)), hx, 'replicate');
gradmagk(:,:,kin) = sqrt(Ix.^2 + Iy.^2);
end

gradmag=sum(gradmagk,3);
    I = gradmag;
     if NR==1
     I=uint16(65535*mat2gray(localnormalize(double(I),75,75)));
     end
    
        Img1=I;
        mean_s=NCHmax/2;  
        mins=0.51;
        maxs=NCHmax/2;
        kernel_s=round(maxs*5);
    
IM4=IO; 
Img=Img1;

        h = fspecial('disk', mean_s);
        Img2 = imfilter(Img,h,'replicate');

        Img3=imcomplement(Img2);
        Img8=imcomplement(Img);

        max3=max(Img3(:));
        min3=min(Img3(:));

        type=16;
     
        types=linspace(mins,maxs,type);
              Chd= uint16(round(type-1)*mat2gray(Img3))+1;
              Imgn=[];
                   for k=1:type
                  
                      h = fspecial('gaussian',[kernel_s kernel_s],types(k));
                                       Imgn(:,:,k) = imfilter(Img8,h,'replicate');
                   end
              adapt_image=zeros(size(Img8));
              
              for k=1:type
                        temp=Imgn(:,:,k);
                      adapt_image(Chd==k)=temp(Chd==k);
                                  
              end

        Img1=uint16(65535*mat2gray(adapt_image));
        
            imwrite(Img1,[nametex 'Adaptive_gaussian_filtered' '.png']);
                 
       Im2=Img1;
%         Im4=uint16(65535*mat2gray(IM4));
       idx=1; 
       tol(idx)=std(double(Im2(:)))/sqrt(3);
       
    MIJ.createImage(Im2);
    MIJ.run('Find Maxima...', sprintf('noise=%d output=[Segmented Particles]',tol(idx)));
    I = MIJ.getCurrentImage;
    
    MIJ.closeAllWindows;
    Img11=uint16(65535*mat2gray(I));
    classmap5=bwlabel(Img11);

    Img12=bwareaopen(Img11>0,MAXS);
    Img13=(Img11>0)-Img12;

lblImg = bwlabel(Img13);

classmap5=bwlabel(Img13);
figure,imshow(label2rgb(classmap5,'jet','k','shuffle'));
classmap=classmap5;

 pcells=unique(classmap);
Vmax=max(IO(:));
         for nk=2:length(pcells)  
              val=pcells(nk);  
            object2 = classmap==val;            
            IGG=mean(IO(object2))
                                          if IGG<(Vmax*LTH)
                                                classmap(object2)=0;

                                          end
         end

LABELM=classmap;
figure,imshow(label2rgb(classmap,'jet','k','shuffle'));
pcells=unique(LABELM);
se=ones(3);
Img13=LABELM>0;
objectbor_map=Img13-imerode(Img13,ones(3));
   LCOLOR5=uint16(cat(3,65535*mat2gray(IO),65535*mat2gray(IO),65535*mat2gray(IO)));                                    
mult=[1 1 1];
                          for ind=1:1
                          col_img2a=mult(ind)*LCOLOR5(:,:,ind);
%                           col_img2a(LABELM==0)=0;  
                          col_img2a(objectbor_map>0)=65535;    
                          LCOLOR5(:,:,ind)=col_img2a;
                          end   

                imwrite(uint16(65535*mat2gray(LCOLOR5)),[nametex 'Segmenation.png']);
                imwrite(uint16(65535*mat2gray(objectbor_map)),[nametex 'Border.png']);


end
      