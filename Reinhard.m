% Getting the Input image
Img = hdrread('memorial.hdr');
[r,c,h] = size(Img);

% Luminance matrix
Lum = zeros(r,c);
for i=1:r
    for j=1:c
        Lum(i,j) = Img(i,j,1)*0.299 + Img(i,j,2)*0.587 + Img(i,j,3)*0.114;
    end
end

% Performing the required operations on the Luminance values
delta = 0.00001;
sum = 0;
for i=1:r
    for j=1:c
        sum = sum + log(delta + Lum(i,j));
    end
end

Lw = exp(sum/(r*c));
L = zeros(r,c);
a = 0.045;

for i=1:r
    for j=1:c
        L(i,j) = (a/Lw)*Lum(i,j);
    end
end

% Defining the constants as directed in the paper
alpha1 = 1/(2*(sqrt(2)));
alpha2 = 1.6*alpha1;
epsilon = 0.05;
phi = 15.0;

R1 = zeros(r,c,9);
R2 = zeros(r,c,9);

for x=1:r
    for y=1:c
        for s=1:9
            sS = 1.6^(s-1);
            R1(x,y,s) = (1/(pi*(alpha1*sS)^2))*(exp(-(x^2 + y^2)/(alpha1 * sS)^2));
            R2(x,y,s) = (1/(pi*(alpha2*sS)^2))*(exp(-(x^2 + y^2)/(alpha2 * sS)^2));
        end
    end
end

% Performing the required Fourier Transforms
for s=1:9
    V1(:,:,s) = ifft(fft(L).*fft(R1(:,:,s)));
    V2(:,:,s) = ifft(fft(L).*fft(R1(:,:,s)));
end

V = zeros(r,c,9);
for x=1:r
    for y=1:c
        for s=1:9
            sS = 1.6^(s-1);
            V(x,y,s) = (V1(x,y,s)-V2(x,y,s))/(((2^phi)*((a)/(sS*sS)))+V1(x,y,s));
        end
    end
end

Ld = zeros(r,c);

for i=1:r
    for j=1:c
        
        sS = 1;
        for s=1:9
            
            if (abs(V(i,j,s)))<epsilon
                break;
            end
            
            if s~=9
                sS = sS * 1.6;
            end
        end
        
        p = 1 + round(log(sS)/log(1.6));
        Ld(i,j) = L(i,j)/(1+V1(i,j,p));
        
    end
end

% Getting the output Image using the transformations done earlier
ImgOut = zeros(r,c,h);
for i = 1:r
    for j = 1:c
        for k =1:h
            ImgOut(i,j,1) = (((Img(i,j,1))/(Lum(i,j))))*(Ld(i,j));
            ImgOut(i,j,2) = (((Img(i,j,2))/(Lum(i,j))))*(Ld(i,j));
            ImgOut(i,j,3) = (((Img(i,j,3))/(Lum(i,j))))*(Ld(i,j));
        end
    end
end
imshow(ImgOut);

% The unscaled Gamma function
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

% The scaled gamma function
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
