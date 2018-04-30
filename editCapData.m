clear

load cap1005FullWithExtra.mat
capInfoBak = capInfo;

capInfo = cell(1,4);
temp=capInfoBak{1};

capInfo{1} = temp(1:82);
capInfo{1}{83}='O9';
capInfo{1}{84}='O10';
% for i=85:396, capInfo{1}{i}=temp(i-2); end
for i=85:396, capInfo{1}{i}=temp{i-2}; end

for k=2:4
    
    capInfo{k} = zeros(396,1);
    capInfo{k}(1:82) = capInfoBak{k}(1:82);
    capInfo{k}(83) = capInfoBak{k}(83);
    capInfo{k}(84) = capInfoBak{k}(85);
    capInfo{k}(85:396) = capInfoBak{k}(83:394);
    
end

save('cap1005FullWithExtra.mat','capInfo');