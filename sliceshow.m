function sliceshow(img,pos,color,clim,label,figName)
% sliceshow(img,pos,color,clim,label,figName)
%
% sliceshow displays a 3D volume allowing the user to click on different slices.
%
% If pos is given as input, then the function will start the display with
% the slices intersecting this position. Defaults to the middle of the
% volume.
%
% color is the colormap the display will use. Defaults to 'jet'.
%
% clim specifies the color limits and defaults to [min max].
%
% If label is given as input, then a colorbar with that label is added.
%
% (c) Nov 02, 2017, Lucas C Parra

if nargin<1 || isempty(img)
    error('At least give us a volume to display')
else
    mydata.img = img;
end

if nargin<2 || isempty(pos)
    mydata.pos = round(size(img)/2);
else
    mydata.pos=pos;
end

if nargin<3 || isempty(color)
    mydata.color = 'jet';
else
    mydata.color = color;
end

if nargin<4 || isempty(clim)
    if min(img(:))==max(img(:))
        mydata.clim = 'auto';
    else
        mydata.clim = [min(img(:)) max(img(:))];
    end
else
    mydata.clim = clim;
end

if nargin<5
    mydata.label=[];
else
    mydata.label=label;
end

if nargin<6
    figName = '';
end

% link the calback function to new figure
fh = figure('WindowButtonDownFcn',@myCallback,'Name',figName,'NumberTitle','off');

% store data as figure property
set(fh,'UserData',mydata);

% first display
showimages


end

function myCallback(fh,~)
% This selects new position depending on which subplot a button is clicked.
% Updates display information in mydata.

mydata = get(fh,'UserData');
pos    = get(gca,'CurrentPoint');
pos    = round(pos(1,1:2));
switch gca
    case mydata.h(1); mydata.pos([1 3]) = pos;
    case mydata.h(2); mydata.pos([2 3]) = pos;
    case mydata.h(3); mydata.pos([1 2]) = pos;
end

% if new position is valid, update mydata as figure property and display
% otherwise do nothing
if ~sum( mydata.pos<1 | mydata.pos>size(mydata.img) )
    set(fh,'UserData',mydata);
    showimages
end

end

function showimages
% Displays the image based on current mydata.

% grab the stuff to display from figure
mydata  = get(gcf,'UserData');

% show the images and keep track of the axes that are generated
h(1)=subplot(2,2,1); imagesc(squeeze(mydata.img(:,mydata.pos(2),:))'); d(1,:)=[1 3];
h(2)=subplot(2,2,2); imagesc(squeeze(mydata.img(mydata.pos(1),:,:))'); d(2,:)=[2 3];
h(3)=subplot(2,2,3); imagesc(squeeze(mydata.img(:,:,mydata.pos(3)))'); d(3,:)=[1 2];

% some aesthetics
val = mydata.img(mydata.pos(1),mydata.pos(2),mydata.pos(3));
for i=1:3
    subplot(2,2,i);
    hold on; plot(mydata.pos(d(i,1)),mydata.pos(d(i,2)),'*','color',[1 1 1]); hold off;
    axis xy; axis equal; axis tight; caxis(mydata.clim);
    title([num2str(mydata.pos(d(i,:))) ': ' num2str(val,'%.2f')])
end

if ~isempty(mydata.label)
    h(4) = subplot(2,2,4); axis off; caxis(mydata.clim);
    h(5) = colorbar('west');
    set(h(5),'YAxisLocation','right','FontSize',18);
    title(h(5), mydata.label,'FontSize',18);
end

colormap(mydata.color);

mydata.h=h;
set(gcf,'UserData',mydata);

end