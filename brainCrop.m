% Crops the image to the brain using segmentation
% Returns the upper and lower bounds of the brain in each direction
% column 1 = R/L, column 2 = A/P, column 3 = S/I
% (c) 2024 Gavin Hsu and Andrew Birnbaum

function bbox = brainCrop(mask)
    wm = mask == 1;
    thresh = [0,0,100];
    [bbox(1,1),bbox(2,1)]=bounds(find(sum(sum(wm,2),3)>thresh(1)));
    [bbox(1,2),bbox(2,2)]=bounds(find(sum(sum(wm,1),3)>thresh(2)));
    [bbox(1,3),bbox(2,3)]=bounds(find(sum(sum(wm,1),2)>thresh(3)));
end