function ridgeLogical = Movie_followRidge(TopographicsImage, startingPixel, scaningRange, GradThresholdMinMax)
%ridgeLogical = Movie_followRidge(TopographicsImage, startingPixel, scaningRange, GradThresholdMinMax)
%
%followRidge will follow the topographic ridge upon which the starting
%pixel is sitting on
%GradThresholdMinMax - in the case of a scalar the function will find the
%minimal gradiant and follow it. In case of it representing the strings
%"min" or "max" it will find to closest minimal or maximal value in the
%following row respectivley.

S = size(TopographicsImage);
ridgeLogical = zeros(S);      % here the ridge will be marked with 1s.

ridgeLogical(1,startingPixel) = 1;        % mark the first point of the the ridge

currentPixel = startingPixel;

if ischar(GradThresholdMinMax)
    %--- methos of following maxima of each row
    
    for rowNum = 2:S(1)
        nearPixelsVector = (currentPixel-scaningRange):(currentPixel+scaningRange);   % this vector indicates the pixels to check for the ridge
        nearPixelsVector(nearPixelsVector>S(2)) = S(2);     % so the indecies won't extend the image
        nearPixelsVector(nearPixelsVector<1) = 1;       	% so the indecies won't extend the image
        
        rowSample = TopographicsImage(rowNum,nearPixelsVector);  	% sample the image in the close2pixeld range
        if strcmpi(GradThresholdMinMax,'min')
            [~, Indx] = min(rowSample);                   % find the closest minima
        elseif strcmpi(GradThresholdMinMax,'max')
            [~, Indx] = max(rowSample);                   % find the closest maxima
        else
            error('GradThresholdMinMax variable incorrect');
        end
        
        %--- if there is noe near change, continue straight:
        if nnz(rowSample-rowSample(1))~=0
            currentPixel = nearPixelsVector(Indx);              % change current pixel for next round
        end
        ridgeLogical(rowNum,currentPixel)=1;
        
    end
    %--- method of following minimal gradient of each row transition
elseif isdouble(GradThresholdMinMax)
    for rowNum = 2:S(1)
        nearPixelsVector = (currentPixel-scaningRange):(currentPixel+scaningRange);   % this vector indicates the pixels to check for the ridge
        nearPixelsVector(nearPixelsVector>S(2)) = S(2);     % so the indecies won't extend the image
        nearPixelsVector(nearPixelsVector<1) = 1;       	% so the indecies won't extend the image
        
        rowSample = TopographicsImage(rowNum,nearPixelsVector);  	% sample the image in the close2pixeld range
        DiffFromPix = rowSample - TopographicsImage(rowNum-1,currentPixel); %vector of gradient in positive time direction
        
        [minGradVal, Indx] = min(abs(DiffFromPix));                   % find the closest maxima
        if minGradVal > GradThresholdMinMax
            break
        end
        
        %--- if there is noe near change, continue straight:
        if nnz(DiffFromPix)~=0
            currentPixel = nearPixelsVector(Indx);              % change current pixel for next round
        end
        ridgeLogical(rowNum,currentPixel)=1;
    end
end

end