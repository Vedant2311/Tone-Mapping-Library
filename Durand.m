% Getting the HDR input file
Img = hdrread('memorial.hdr');
[r,c,h] = size(Img);

% Luminance matrix
Lum = zeros(r,c);
for i=1:r
    for j=1:c
        Lum(i,j) = Img(i,j,1)*0.299 + Img(i,j,2)*0.587 + Img(i,j,3)*0.114;
    end
end

% Computing the Gradient matrix as directed in the paper
H = zeros(r,c);
for i=1:r
    for j = 1:c
        H(i,j) = log10(Lum(i,j));
    end
end

gradH = zeros(r,c,2);
for i=1:r
    for j=1:c
        if i==1 
           A = H(i+1,j);
           B=0;        
        elseif i==r
            A = 0;
            B= H(i-1,j);
        else
            A= H(i+1,j);
            B = H(i-1,j);
        end
   
        if j==1   
            C = H(i,j+1);
            D = 0;
        elseif j==c
            C = 0;
            D = H(i,j-1);
        else
            C = H(i,j+1);
            D = H(i,j-1);
        end
        
        x = (A-B)/2;
        y = (C-D)/2;
        gradH(i,j,:) = [x,y];
       
    end
end

sum = 0;
for i=i:r
    for j=1:c   
        sum = sum + sqrt(gradH(i,j,1)^2 + gradH(i,j,2)^2);
    end
end

% Getting the constants as directed in the paper
average = sum/r*c;
alpha = 0.1*average;
beta = 0.85;

% Performing the transformations as directed in the paper
phi = zeros(r,c);
for i=1:r
    for j=1:c
        phi(i,j) = (alpha/sqrt(gradH(i,j,1)^2 + gradH(i,j,2)^2)) * (sqrt(gradH(i,j,1)^2 + gradH(i,j,2)^2)/alpha)^(beta);
    end
end

G = gradH .* phi;
divG = zeros(r,c);
for i=1:r
    for j=1:c
        
        if ((i==1) && (j==1))
            
            divG(i,j) = G(i,j,1) + G(i,j,2);
        elseif i==1
            
            divG(i,j) = G(i,j,1) + G(i,j,2) - G(i,j-1,2);
        elseif j==1
            
            divG(i,j) = G(i,j,1) - G(i-1,j,1) + G(i,j,2);
        else
            
            divG(i,j) = G(i,j,1) - G(i-1,j,1) + G(i,j,2) - G(i,j-1,2);
        end
    end
end
divG = fillmissing(divG,'previous');
f = fft(divG);

fi = zeros(r,c);
Lx = r;
Ly = c;
for i=1:r
    for j = 1:c
        fi(i,j) = (f(i,j))/(((2*pi*(i/Lx)*(1j))^2)+((2*pi*(j/Ly)*(1j))^2));
    end
end
I = real(ifft(fi));

for i=1:r
    for j=1:c
        I(i,j) = 10^(I(i,j));
        if j==1 
            I(i,j) = 0;
        end 
    end
end

% Getting the output image after all the transformations done earlier
ImgOutLum = zeros(r,c,h);
for i = 1:r
    for j = 1:c
        for k =1:h
            ImgOutLum(i,j,1) = (((Img(i,j,1))/(Lum(i,j))))^(0.5)*(I(i,j));
            ImgOutLum(i,j,2) = (((Img(i,j,2))/(Lum(i,j))))^(0.5)*(I(i,j));
            ImgOutLum(i,j,3) = (((Img(i,j,3))/(Lum(i,j))))^(0.5)*(I(i,j));
        end
    end
end

imshow(ImgOutLum);
