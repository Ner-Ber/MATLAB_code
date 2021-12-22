function h = pyrHeight(im, maxLevels)
%pyrHeight calculates the wanted pyramid height
% h - will be the pyramid height
% s - vector for image size

    s = size(im);
    n = floor(log2((min(s)/16)));           % calculates number ofiteration needed to reach 16 pix
    h = min(maxLevels, n);                  % wanted height of pyramid

end