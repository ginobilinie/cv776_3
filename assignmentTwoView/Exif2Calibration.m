function [K]=Exif2Calibration(exif)
%
%
%

if isfield(exif,'FocalLengthIn35mmFilm')
    
        focal=exif.CPixelXDimension*exif.FocalLengthIn35mmFilm/35;
        K = [focal 0 exif.CPixelXDimension/2;
             0 focal exif.CPixelYDimension/2;
             0 0 1];
        return;
else
        K=eye(3)
        display('Unknown EXIF information');
end
        