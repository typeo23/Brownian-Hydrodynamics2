struveMat = zeros(1,length(1:1e-2:100));
struveMex = struveMat;

ind = 1;
tic
for i = 1e-5:1e-2:1000
    struveMat(ind) = bessely(0,i);
    ind = ind+1;
end
toc
ind = 1;
tic
for i = 1e-5:1e-2:1000
    struveMex(ind) = besselyn(0,i);
    ind = ind+1;
end
toc