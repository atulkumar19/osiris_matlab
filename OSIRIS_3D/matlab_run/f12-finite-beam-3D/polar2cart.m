function [x,y]= polar2cart (mag, ang_in_deg)
x = mag * cos(ang_in_deg*pi/180);
y = j * mag * sin(ang_in_deg*pi/180);
 

function [r, ar, ad] = cart2polar(x)
r = abs(x);
ar = angle(x);
ad = ar*180/pi;

And now we test them: 

% Clear memory and screen. Avoid double-blank lines
clear; clc; format compact 

[x, y] = polar2cart(2, 30.5)
[r, ar, ad] = cart2polar(7 + 18i)
[r, ar, ad] = cart2polar(0 - 46.8i)
 