# Tone-Mapping-Library

Real-world scenes often have a very large dynamic range (the ratio of brightest to darkest intensities) which can span several orders of magnitude. Such high dynamic range (HDR) images cannot be reproduced directly on conventional displays. For a more naturalistic appearance, the range of intensities has to be compressed to the low dynamic range of the display, while approximately maintaining the appearance of the image. This process is known as tone mapping or dynamic range compression.

## Linear and Logarithmic rescaling

Loaded an HDR image into the program and visualized it by linearly rescaling the pixel values. Different ranges for linearly rescaling can be done and not just **0-255**

Note that the HDR image contains linear intensity values, while low dynamic range images are quantized nonlinearly, so applied a gamma correction to get a reasonable image. Created two kinds of Gamma functions: 

  1. Unscaled: Applied the function without scaling the input image to 0-1
  2. Scaled: Applied the function after the scaling, and then applying linear rescaling to **0-255** or **max(0,minOrig)-min(255,maxOrig)**

As a baseline tone mapping algorithm, perform rescaling in the log-luminance domain. That is, compute the luminance L=0.299R+0.587G+0.114B and take its log, rescale it to decrease the dynamic range, then undo the log to recover luminances, and use the original colour ratios R/L,G/L,B/L to reconstruct a colour image. Choose a scaling in the log domain to get an output dynamic range of about 100:1.

