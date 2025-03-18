function alignImagesYoPro_BloodVesssels(yoproIm, bvImage)

g = read_Tiffs(yoproIm);
r = read_Tiffs(bvImage);
yoproImSize = size(g);
bvImageSize = size(r);

diffSize = yoproImSize- bvImageSize;

if diffSize(1)>0
    r = paddata(r, [size(r,1)+diffSize(1) 0] , Side="leading");
else
    g = paddata(g, [size(g,1)+abs(diffSize(1)) 0] , Side="leading");
end

if diffSize(2)>0
    r = paddata(r, [0 size(r,2)+diffSize(2)] , Side="leading");
else
    g = paddata(g, [0 size(g,2)+abs(diffSize(2))] , Side="leading");
end

imFuse = cat(3,r,g, zeros(size(g)));
savePath = regexprep(yoproIm, '_C[12]_', '_');
options.color = true;
saveastiff(imFuse,savePath,options);
end