function stdArray = stdGPU(dat, dim)

if nargin < 2 || isempty(dim)
    dim =  ndims(dat);
end

datSz = size(dat);

% get reshape inputs
reshapeTxt = ['[],' num2str(size(dat,dim)) ];

% reshape data for ease
eval(['dat = reshape(dat,' reshapeTxt ');']);

%std calculation
mean_x = uint16(sum(dat,2)/size(dat,2));
xc = dat - mean_x;
stdArray = sqrt(sum(xc .* xc, 2)) / sqrt(size(dat,2) - 1);

reshapeIndx = datSz;
reshapeIndx(dim) = [];

stdArray = reshape(stdArray,reshapeIndx);

end