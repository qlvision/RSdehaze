function newim= rsshow(band1,band2,band3,perc)
% show rs image using linear stretch, ENVI style
if nargin<4
    I1 = uint8( 255*fnnorm( double(band1) ) );
    I2 = uint8( 255*fnnorm( double(band2) ) );
    I3 = uint8( 255*fnnorm( double(band3) ) );
    newim = cat(3,I3,I2,I1);
    imshow( newim );
else
    b1dist = rebound( double(band1) ,perc);
    b2dist = rebound( double(band2) ,perc);
    b3dist = rebound( double(band3) ,perc);
    I1 = uint8( 255*fnnorm(b1dist) );
    I2 = uint8( 255*fnnorm(b2dist) );
    I3 = uint8( 255*fnnorm(b3dist) );
    newim = cat(3,I3,I2,I1);
    imshow( newim );
end

end

function bd = rebound(band,perc)
    [y_size,x_size] = size(band);
    nhist = y_size*x_size;
    bmin = min(band(:));
    bmax = max(band(:));
    x= bmin:0.5:bmax;
    bhist = hist(band(:),x);
    sum1 = 0;
    for i=1:length(x)
        sum1 = sum1 + bhist(i);
        if sum1/nhist > perc/100
            blow = x(i);
            break;
        end
    end
    sum1 = 0;
    for i=length(x):-1:1
        sum1 = sum1 + bhist(i);
        if sum1/nhist > perc/100
            bhig = x(i);
            break;
        end
    end
    bd = min(max(band,blow),bhig);
end

function imnorm = fnnorm(im,bmin,bmax)
im = double(im);
if nargin==1
bmin=min(min(min(im)));
bmax=max(max(max(im)));
end
bmin = double(bmin);
bmax = double(bmax);
eps = 0.000001;
imnorm = (im-bmin+eps)./(bmax-bmin+eps);
end