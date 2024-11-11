function vec = str2vec(inputstr)

M = size(inputstr, 1);
vec = [];
	
for ii =1:M
    vec = [vec; str2num(deblank(inputstr(ii,:)))];
end
