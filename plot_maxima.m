function [ImgOut]=plot_maxima(Imgin,points,se)                         

edgeflag3 = imdilate(points,se);
                            
                            a1=Imgin(:,:,1);
                            a1(edgeflag3>0)=2^16-1;
                            
                            a2=Imgin(:,:,1);
                            a2(edgeflag3>0)=0;
                            
                            a3=Imgin(:,:,1);
                            a3(edgeflag3>0)=0;
                            
                            ImgOut(:,:,1)=a1;
                            ImgOut(:,:,2)=a2;
                            ImgOut(:,:,3)=a3;