Img = hdrread('memorial.hdr');
[r,c,h] = size(Img);
maxRGB = max(max(Img));
minRGB = min(min(Img));

ImgOutSimple = zeros(r,c,h);
ImgOutLow = zeros(r,c,h);
ImgOutHigh = zeros(r,c,h);
ImgOutGamma = zeros(r,c,h);

% Linear rescaling

for i=1:r
    for j = 1:c
        for k = 1:h
            ImgOutSimple(i,j,k) = ((Img(i,j,k)- minRGB(k))/(maxRGB(k)-minRGB(k)))*255;
            ImgOutLow(i,j,k) = ((Img(i,j,k)- minRGB(k))/(maxRGB(k)-minRGB(k)))*100;
            ImgOutHigh(i,j,k) = ((Img(i,j,k)- minRGB(k))/(maxRGB(k)-minRGB(k)))*1000;
            
        end
    end
end

% Writing the Linearly rescaled images

%imwrite(uint8((ImgOutSimple)),'C:\Users\Lenovo\Downloads\COL783_Assgn\A1\memorialSimple.jpg','jpg');
%imwrite(uint8((ImgOutLow)),'C:\Users\Lenovo\Downloads\COL783_Assgn\A1\memorialSimpleLow.jpg','jpg');
%imwrite(uint8((ImgOutLow)),'C:\Users\Lenovo\Downloads\COL783_Assgn\A1\memorialSimpleHigh.jpg','jpg');

% Applying the Gamma correction for the Linearly rescaled value of the
% range 0 - 255

for i=1:r
    for j = 1:c
        for k = 1:h
           if (ImgOutSimple(i,j,k)/255)<=0.0031308 
              ImgOutGamma(i,j,k) = 12.92*(ImgOutSimple(i,j,k)/255); 
           else
               ImgOutGamma(i,j,k) =  (((ImgOutSimple(i,j,k)/255)^(1/2.4))*1.055 + (-0.055))*1;
           end
        
           ImgOutGamma(i,j,k) = ImgOutGamma(i,j,k)*255;
           
        end
    end
end

% Writing the gamma correction output
%imwrite(uint8((ImgOutGamma)),'C:\Users\Lenovo\Downloads\COL783_Assgn\A1\memorialGamma_1.jpg','jpg');

% Luminance matrix 
ImgOutLum = zeros(r,c,h);
Lum = zeros(r,c);

for i=1:r
    for j=1:c
        Lum(i,j) = Img(i,j,1)*0.299 + Img(i,j,2)*0.587 + Img(i,j,3)*0.114;
    end
end

% Log luminance matrix
Lumlog = zeros(r,c);

for i=1:r
    for j = 1: c
        Lumlog(i,j) = log10(Lum(i,j));
    end
end

maxVal = max(max(Lumlog));
minVal = min(min(Lumlog));

% Scaling the Log luminance to the range 0 - 2 => Luminance is scaled as
% 1 - 100

for i=1:r
    for j = 1: c
        Lumlog(i,j) = ((Lumlog(i,j)-minVal)/(maxVal-minVal))*2-1;
        Lumlog(i,j) = 10^(Lumlog(i,j));
    end
end

% Using the constancy of R/L, B/L, G/L to recover the colored image from
% the scaled Luminance matrix

for i = 1:r
    for j = 1:c
        for k =1:h
            ImgOutLum(i,j,1) = (((Img(i,j,1))/(Lum(i,j))))*(Lumlog(i,j));
            ImgOutLum(i,j,2) = (((Img(i,j,2))/(Lum(i,j))))*(Lumlog(i,j));
            ImgOutLum(i,j,3) = (((Img(i,j,3))/(Lum(i,j))))*(Lumlog(i,j));
        end
    end
end

% Writing the corrected image
imshow(ImgOutLum);
%imwrite(uint8((ImgOutLum)),'C:\Users\Lenovo\Downloads\COL783_Assgn\A1\memorial_LogLumi.jpg','jpg');
     
%Convolution matrices...

A = 1.5;
Sharpen = [0,-1,0;-1,5,-1;0,-1,0];
Laplace = [0,1,0; 1 , -4,1; 0,1,0];
Highboost = [-1,-1,-1;-1,A+8,-1;-1,-1,-1];
Sobel = [-1,-2,-1;0,0,0;1,2,1];

% Using the Unsharp Masking technique 
ImgSharpen = convol(ImgOutLum,r,c,Sharpen);
%imshow(unscaledGamma(ImgSharpen,r-2,c-2,h));

% Using the Laplace Technique
ImgLap = convol(ImgOutLum, r,c,Laplace);
%imshow(unscaledGamma(ImgLap,r-2,c-2,h));

% Using the Highboost filtering
ImgBoost = convol(ImgOutSimple,r,c,Highboost);
%imshow(ImgBoost);


% Using the Gradient filter
ImgGrad = convol(ImgOutLum,r,c,Sobel);
%imshow(unscaledGamma(ImgGrad,r-2,c-2,h));

% Using contrast stretching
%imshow(contrast(ImgOutLum,r,c,h,255/3,2*255/3,255/6,5*255/6));
%imshow(unscaledGamma(contrast(ImgOutSimple,r,c,h,255/3,2*255/3,255/6,5*255/6),r,c,h));

%Using Histogram equalization
%imshow(unscaledGamma(histo(ImgOutLum,r,c,h),r,c,h));

function A = convol(Bmat,r,c,Cmat)

  Res = zeros(r-2,c-2,3);
  
  for i=1:(r-2)
      for j = 1:(c-2)
         
          for k=1:3
          sum = 0;
         
             for u=0:2
                 for v=0:2
                     sum = sum + Bmat(i+u,j+v,k)*Cmat(3-u,3-v);
                 end
             end

             Res(i,j,k) = sum;
             
          end
      end
  end
  
  A = Res;
  
end

function A = gamma(ImgOutSimple,r,c,h)

    ImgOutGamma = zeros(r,c,h);
    
    for i=1:r
        for j = 1:c
            for k = 1:h
               if (ImgOutSimple(i,j,k)/255)<=0.0031308 
                  ImgOutGamma(i,j,k) = 12.92*(ImgOutSimple(i,j,k)/255); 
               else
                   ImgOutGamma(i,j,k) =  (((ImgOutSimple(i,j,k)/255)^(1/2.4))*1.055 + (-0.055))*1;
               end

               ImgOutGamma(i,j,k) = ImgOutGamma(i,j,k)*255;

            end
        end
    end

    A = ImgOutGamma;

end

function A = contrast(Img,r,c,h,r1,r2,s1,s2)

  A = zeros(r,c,h);
  
  for i=1:r
      for j=1:c
          for k=1:h
              
              if Img(i,j,k)<r1
                  
                  A(i,j,k)=Img(i,j,k)/s1;
                  
              elseif Img(i,j,k)<r2
                  
                  A(i,j,k) = s1 + ((Img(i,j,k)-r1)/(r2-r1))*(s2-s1);
              else
              
                  A(i,j,k) = s2 + ((Img(i,j,k)-r2)/(255-r2))*(255-s2);
              end
          end
      end
  end
    
end

function A = histo(Img,r,c,h)

    minVal = floor(min(min(Img)));
    maxVal = floor(max(max(Img)));
    arr1 = zeros(maxVal(1)-minVal(1)+1);
    arr2 = zeros(maxVal(2)-minVal(2)+1);
    arr3 = zeros(maxVal(3)-minVal(3)+1);
    
    A = zeros(r,c,h);
    
    for i=1:r
        for j=1:c
            
            arr1(floor(Img(i,j,1))-minVal(1)+1) = arr1(floor(Img(i,j,1))-minVal(1)+1)+1;
            arr2(floor(Img(i,j,2))-minVal(2)+1) = arr2(floor(Img(i,j,2))-minVal(2)+1)+1;
            arr3(floor(Img(i,j,3))-minVal(3)+1) = arr3(floor(Img(i,j,3))-minVal(3)+1)+1;
        end
    end
    
    for i=1:(maxVal(1)-minVal(1)+1)
        
        if i~=1
            arr1(i) = arr1(i-1)+arr1(i);
        end
    end
    
    for i=1:(maxVal(2)-minVal(2)+1)
        
        if i~=1
            arr2(i) = arr2(i-1)+arr2(i);
        end
    end
    
    for i=1:(maxVal(3)-minVal(3)+1)
        
        if i~=1
            arr3(i) = arr3(i-1)+arr3(i);
        end
    end
    
    for i=1:r
        for j=1:c
            
            A(i,j,1) = ((arr1(floor(Img(i,j,1))-minVal(1)+1)-arr1(1))/((r*c)-arr1(1)))*255;
            A(i,j,2) = ((arr2(floor(Img(i,j,2))-minVal(2)+1)-arr2(1))/((r*c)-arr2(1)))*255;
            A(i,j,3) = ((arr3(floor(Img(i,j,3))-minVal(3)+1)-arr3(1))/((r*c)-arr3(1)))*255;
            
        end
    end
    
    
end

function A = unscaledGamma(ImgOutSimple,r,c,h)

    ImgOutGamma = zeros(r,c,h);
    
    for i=1:r
        for j = 1:c
            for k = 1:h
               if (ImgOutSimple(i,j,k))<=0.0031308 
                  ImgOutGamma(i,j,k) = 12.92*(ImgOutSimple(i,j,k)); 
               else
                   ImgOutGamma(i,j,k) =  (((ImgOutSimple(i,j,k))^(1/2.4))*1.055 + (-0.055))*1;
               end
            end
        end
    end

    A = ImgOutGamma;

end