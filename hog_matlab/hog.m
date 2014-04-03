function H = getHistogram(magnitudes, angles, numBins)
% GETHISTOGRAM Computes a histogram for the supplied gradient vectors.
%   H = getHistogram(magnitudes, angles, numBins)
%   
%   This function takes the supplied gradient vectors and places them into a
%   histogram with 'numBins' based on their unsigned orientation. 
%
%   "Unsigned" orientation means that, for example, a vector with angle 
%   -3/4 * pi will be treated the same as a vector with angle 1/4 * pi. 
%
%   Each gradient vector's contribution is split between the two nearest bins,
%   in proportion to the distance between the two nearest bin centers.
%   
%   A gradient's contribution to the histogram is equal to its magnitude;
%   the magnitude is divided between the two nearest bin centers.
%
%   Parameters:
%     magnitudes - A column vector storing the magnitudes of the gradient 
%                  vectors.
%     angles     - A column vector storing the angles in radians of the 
%                  gradient vectors (ranging from -pi to pi)
%     numBins    - The number of bins to place the gradients into.
%   Returns:
%     A row vector of length 'numBins' containing the histogram.

% $Author: ChrisMcCormick $    $Date: 2013/11/27 22:00:00 $    $Revision: 1.1 $

% Revision Notes:
%  v1.1
%    - The function was actually hardcoded to 9 bins; it now properly supports
%      specifying 'numBins'.
%    - It now returns an unsigned histogram. This has been shown to improve 
%      accuracy for person detection.

% Compute the bin size in radians. 180 degress = pi.
binSize = pi / numBins;

% The angle values will range from 0 to pi.
minAngle = 0;

% Make the angles unsigned by adding pi (180 degrees) to all negative angles.
angles(angles < 0) += pi;

% The gradient angle for each pixel will fall between two bin centers.
% For each pixel, we split the bin contributions between the bin to the left
% and the bin to the right based on how far the angle is from the bin centers.

% For each pixel's gradient vector, determine the indeces of the bins to the
% left and right of the vector's angle.
%
% The histogram needs to wrap around at the edges--vectors on the far edges of
% the histogram (i.e., close to -pi or pi) will contribute partly to the bin
% at that edge, and partly to the bin on the other end of the histogram.
% For vectors with an orientation close to 0 radians, leftBinIndex will be 0. 
% Likewise, for vectors with an orientation close to pi radians, rightBinIndex
% will be numBins + 1. We will fix these indeces after we calculate the bin
% contribution amounts.
leftBinIndex = round((angles - minAngle) / binSize);
rightBinIndex = leftBinIndex + 1;

% For each pixel, compute the center of the bin to the left.
leftBinCenter = ((leftBinIndex - 0.5) * binSize) - minAngle;

% For each pixel, compute the fraction of the magnitude
% to contribute to each bin.
rightPortions = angles - leftBinCenter;
leftPortions = binSize - rightPortions;
rightPortions = rightPortions / binSize;
leftPortions = leftPortions / binSize;

% Before using the bin indeces, we need to fix the '0' and '10' values.
% Recall that the histogram needs to wrap around at the edges--bin "0" 
% contributions, for example, really belong in bin 9.
% Replace index 0 with 9 and index 10 with 1.
leftBinIndex(leftBinIndex == 0) = numBins;
rightBinIndex(rightBinIndex == (numBins + 1)) = 1;

% Create an empty row vector for the histogram.
H = zeros(1, numBins);

% For each bin index...
for (i = 1:numBins)
    % Find the pixels with left bin == i
    pixels = (leftBinIndex == i);
        
    % For each of the selected pixels, add the gradient magnitude to bin 'i',
    % weighted by the 'leftPortion' for that pixel.
    H(1, i) += sum(leftPortions(pixels)' * magnitudes(pixels));
    
    % Find the pixels with right bin == i
    pixels = (rightBinIndex == i);
        
    % For each of the selected pixels, add the gradient magnitude to bin 'i',
    % weighted by the 'rightPortion' for that pixel.
    H(1, i) += sum(rightPortions(pixels)' * magnitudes(pixels));
end    

end


function H = getHOGDescriptor(img)
% GETHOGDESCRIPTOR computes a HOG descriptor vector for the supplied image.
%   H = getHOGDescriptor(img)
%  
%   This function takes a 130 x 66 pixel gray scale image (128 x 64 with a 
%   1-pixel border for computing the gradients at the edges) and computes a 
%   HOG descriptor for the image, returning a 3,780 value column vector.
%
%   Parameters:
%     img - A grayscale image matrix with 130 rows and 66 columns.
%   Returns:
%     A column vector of length 3,780 containing the HOG descriptor.
%
%    The intent of this function is to implement the same design choices as
%    the original HOG descriptor for human detection by Dalal and Triggs.
%    Specifically, I'm using the following parameter choices:
%      - 8x8 pixel cells
%      - Block size of 2x2 cells
%      - 50% overlap between blocks
%      - 9-bin histogram
%
%     The above parameters give a final descriptor size of 
%     7 blocks across x 15 blocks high x 4 cells per block x 9 bins per hist
%     = 3,780 values in the final vector.
%
%    A couple other important design decisions:
%    - Each gradient vector splits its contribution proportionally between the
%      two nearest bins
%    - For the block normalization, I'm using L2 normalization.
%
%    Differences with OpenCV implementation:
%    - OpenCV uses L2 hysteresis for the block normalization.
%    - OpenCV weights each pixel in a block with a gaussian distribution
%      before normalizing the block.
%    - The sequence of values produced by OpenCV does not match the order
%      of the values produced by this code.

% $Author: ChrisMcCormick $    $Date: 2013/11/27 22:00:00 $    $Revision: 1.1 $

% Revision Notes:
%  v1.1
%    - Changed the 'hy' filter mask to match the OpenCV filter mask (it's now
%      [-1; 0; 1] instead of [1; 0; -1]).
%    - Added a required one-pixel border around the image for computing
%      gradients for the edge pixels.

% The number of bins to use in the histograms.
numBins = 9;

% The cells are 8 x 8 pixels.
cellSize = 8;

% Empty vector to store computed descriptor.
H = [];

% Verify the image size is 66 x 130.
[height, width] = size(img); 
if ((width ~= 66) || (height ~= 130))
    disp("Image size must be 130 x 66 pixels (128x64 with 1px border).\n");
    return;
end

% Compute the number cells horizontally and vertically (should be 8 x 16).    
numHorizCells = 8;
numVertCells = 16;

% ===============================
%    Compute Gradient Vectors
% ===============================
% Compute the gradient vector at every pixel in the image.

% Create the operators for computing image derivative at every pixel.
hx = [-1,0,1];
hy = hx';

% Compute the derivative in the x and y direction for every pixel.
dx = imfilter(double(img), hx);
dy = imfilter(double(img), hy);

% Remove the 1 pixel border.
dx = dx(2 : (rows(dx) - 1), 2 : (columns(dx) - 1));
dy = dy(2 : (rows(dy) - 1), 2 : (columns(dy) - 1)); 

% Convert the gradient vectors to polar coordinates (angle and magnitude).
angles = atan2(dy, dx);
magnit = ((dy.^2) + (dx.^2)).^.5;

% =================================
%     Compute Cell Histograms 
% =================================
% Compute the histogram for every cell in the image. We'll combine the cells
% into blocks and normalize them later.

% Create a three dimensional matrix to hold the histogram for each cell.
histograms = zeros(numVertCells, numHorizCells, numBins);

% Cast the image to floating point values.
img = double(img);

% For each cell in the y-direction...
for row = 0:(numVertCells - 1)
    
    % Compute the row number in the 'img' matrix corresponding to the top
    % of the cells in this row. Add 1 since the matrices are indexed from 1.
    rowOffset = (row * cellSize) + 1;
    
    % For each cell in the x-direction...
    for col = 0:(numHorizCells - 1)
    
        % Select the pixels for this cell.
        
        % Compute column number in the 'img' matrix corresponding to the left
        % of the current cell. Add 1 since the matrices are indexed from 1.
        colOffset = (col * cellSize) + 1;
        
        % Compute the indices of the pixels within this cell.
        rows = rowOffset : (rowOffset + cellSize - 1);
        cols = colOffset : (colOffset + cellSize - 1);
        
        % Select the angles and magnitudes for the pixels in this cell.
        cellAngles = angles(rows, cols); 
        cellMagnitudes = magnit(rows, cols);
    
        % Compute the histogram for this cell.
        % Convert the cells to column vectors before passing them in.
        histograms(row + 1, col + 1, :) = getUnsHistogram(cellMagnitudes(:), cellAngles(:), numBins);
    end
    
end

% ===================================
%       Block Normalization
% ===================================    

% Take 2 x 2 blocks of cells and normalize the histograms within the block.
% Normalization provides some invariance to changes in contrast, which can
% be thought of as multiplying every pixel in the block by some coefficient.

% For each cell in the y-direction...
for row = 1:(numVertCells - 1)    
    % For each cell in the x-direction...
    for col = 1:(numHorizCells - 1)
    
        % Get the histograms for the cells in this block.
        blockHists = histograms(row : row + 1, col : col + 1, :);
        
        % Put all the histogram values into a single vector (nevermind the 
        % order), and compute the magnitude.
        % Add a small amount to the magnitude to ensure that it's never 0.
        magnitude = norm(blockHists(:)) + 0.01;
    
        % Divide all of the histogram values by the magnitude to normalize 
        % them.
        normalized = blockHists / magnitude;
        
        % Append the normalized histograms to our descriptor vector.
        H = [H; normalized(:)];
    end
end
    
end


A=imread("crop_000606b.png");
H=getHOGDescriptor((A(:,:,1)));
imshow(H)
