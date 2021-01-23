# Tone-Mapping-Library

Real-world scenes often have a very large dynamic range (the ratio of brightest to darkest intensities) which can span several orders of magnitude. Such high dynamic range (HDR) images cannot be reproduced directly on conventional displays. For a more naturalistic appearance, the range of intensities has to be compressed to the low dynamic range of the display, while approximately maintaining the appearance of the image. This process is known as tone mapping or dynamic range compression. A number of algorthms which help in Tone mapping have been considered and implemented here in MATLAB. The user is free to choose any of these depending on the requirements of speed, details etc

## Linear and Logarithmic rescaling

Loaded an HDR image into the program and visualized it by linearly rescaling the pixel values. Different ranges for linearly rescaling can be done and not just **0-255**

Note that the HDR image contains linear intensity values, while low dynamic range images are quantized nonlinearly, so applied a gamma correction to get a reasonable image. Created two kinds of Gamma functions: 

  1. Unscaled: Applied the function without scaling the input image to 0-1
  2. Scaled: Applied the function after the scaling, and then applying linear rescaling to **0-255** or **max(0,minOrig) - min(255,maxOrig)**

As a baseline tone mapping algorithm, performed rescaling in the log-luminance domain. That is, compute the luminance **L=0.299R+0.587G+0.114B** and take its log, rescale it to decrease the dynamic range, then undo the log to recover luminances, and use the original colour ratios R/L,G/L,B/L to reconstruct a colour image. Chose a scaling in the log domain to get an output dynamic range of about 100. (Found 0.1 - 10 to be showing the best results in general)

## Detail Enhancement

You should find that details in the image become weaker, because logarithmic rescaling indiscriminately compresses both large-scale intensity variations and local contrast. To counteract this effect, the following techniques are implemented for tone-maping. You can find them in the code **tonemapping.m**:

  1. Sharpening
  2. Laplacian 
  3. Sobel filter
  4. Highboosting
  5. Contrast stretching
  6. Histogram equalization
  
## More involved Tone mapping algorithms

1. Reinhard et al. (http://www.cs.utah.edu/~reinhard/cdrom/tonemap.pdf)
2. Durand et al. (https://people.csail.mit.edu/fredo/PUBLI/Siggraph2002/DurandBilateral.pdf)
