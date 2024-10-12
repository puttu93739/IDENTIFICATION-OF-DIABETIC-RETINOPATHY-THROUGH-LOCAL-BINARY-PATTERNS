function [pixval] = biinterp(gpx,gpy,image)

[rsize,csize] = size(image);
if (gpx<=0) || (gpy<=0) || (gpx>rsize) || (gpy>csize)
    pixval = 0;
    %fprintf('\npixval=%d',pixval);
else
    x1=floor(gpx); x2=ceil(gpx); y1=floor(gpy); y2=ceil(gpy);
    %fprintf('\nx1=%d \tx2=%d \ty1=%d \ty2=%d',x1,x2,y1,y2);
    q11=pixlook(x1,y1,image); q12=pixlook(x1,y2,image); q21=pixlook(x2,y1,image); q22=pixlook(x2,y2,image);
   % fprintf('\nq11=%f \tq12=%f \tq21=%f \tq22=%f',q11,q12,q21,q22);
    fxy1 = ( (((x2-gpx)/(x2-x1))* q11) + (((gpx-x1)/(x2-x1))* q21) );
    fxy2 = ( (((x2-gpx)/(x2-x1))* q12) + (((gpx-x1)/(x2-x1))* q22) );
    %fprintf('\nfxy1=%f \t fxy2=%f',fxy1,fxy2);
    fxy= ( (((y2-gpy)/(y2-y1))* fxy1) + (((gpy-y1)/(y2-y1))* fxy2) );
    pixval = fxy;
    %fprintf('\npixval=%f',pixval);
end
end

function [pixval] = pixlook(x,y,image)
[rsize,csize] = size(image);
if (x<=0) || (y<=0) || (x>rsize) || (y>csize)
    pixval = 0;
else
    pixval=image(x,y);
    %xc=x;yc=y;
end
end