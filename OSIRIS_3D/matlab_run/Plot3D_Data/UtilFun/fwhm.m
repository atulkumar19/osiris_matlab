function width = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)


y = y / max(y);
N = length(y);
lev50 = 0.5;
if y(1) < lev50                  % find index of center (max or min) of pulse
    [~,centerindex]=max(y);
else
    [~,centerindex]=min(y);
end
ii= 2;
while sign(y(ii)-lev50) == sign(y(ii-1)-lev50)
    ii= ii+1;
end                                   %first crossing is between v(ii-1) & v(ii)
interp = (lev50-y(ii-1)) / (y(ii)-y(ii-1));
tlead = x(ii-1) + interp*(x(ii)-x(ii-1));
ii= centerindex+1;                    %start search for next crossing at center
while ((sign(y(ii)-lev50) == sign(y(ii-1)-lev50)) && (ii<= N-1))
    ii= ii+1;
end
if ii~= N
    interp = (lev50-y(ii-1)) / (y(ii)-y(ii-1));
    ttrail = x(ii-1) + interp*(x(ii)-x(ii-1));
    width = ttrail - tlead;
else
    width = NaN;
end
