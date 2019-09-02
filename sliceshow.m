function sliceshow(img,pos,color,clim,label,figName,vecImg,mri2mni)
% sliceshow(img,pos,color,clim,label,figName,vecImg,mri2mni)
%
% sliceshow displays a 3D volume allowing the user to click on different
% slices. If the variable "vecImg" is specified, sliceshow will also
% overlay a 3D vectorial field represented by arrows on the slices.
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
% figName gives a name to the figure window.
% 
% vecImg is the 3D vectorial field, e.g. an electric field.
% 
% mri2mni is the mapping from MRI voxel space to the MNI space, so
% sliceshow can display both the voxel and MNI coordinates.
%
% (c) Nov 02, 2017, Lucas C Parra
% (c) June 29, 2018, Yu (Andy) Huang
% (c) August 2019, Yu (Andy) Huang

if nargin<1 || isempty(img)
    error('At least give us a volume to display')
else
    [Nx,Ny,Nz] = size(img);
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
    if all(isnan(img(:)))
        error('The image volume you provided does not have any meaningful values.');
    elseif double(min(img(:)))==double(max(img(:)))
        mydata.clim = 'auto';
    else
        mydata.clim = double([min(img(:)) max(img(:))]);
    end
else
    mydata.clim = clim;
end

if nargin<5 || isempty(label)
    mydata.label=[];
else
    mydata.label=label;
end

if nargin<6 || isempty(figName)
    figName = '';
end

if nargin<7 || isempty(vecImg)
    mydata.vecImg = [];
else
    temp = size(vecImg);
    if any(temp(1:3)~=[Nx,Ny,Nz]) || temp(4)~=3
        error('Vector field does not have correct size.');
    end
    mydata.vecImg = vecImg;
    [xi,yi,zi] = ndgrid(1:Nx,1:Ny,1:Nz);
    mydata.xi = xi; mydata.yi = yi; mydata.zi = zi;
end

if nargin<8 || isempty(mri2mni)
    mydata.mri2mni = [];
else
    if any(size(mri2mni)~=[4 4]) || any(mri2mni(4,:)~=[0 0 0 1])
        error('Unrecognized format of the voxel-to-MNI mapping.');
    end
    mydata.mri2mni = mri2mni;
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
h(1)=subplot(2,2,1);
imagesc(squeeze(mydata.img(:,mydata.pos(2),:))'); d(1,:)=[1 3];
if ~isempty(mydata.vecImg)
    hold on;
    rngx=1:5:size(mydata.xi,1); rngy=mydata.pos(2); rngz=1:5:size(mydata.xi,3);
    quiver3(mydata.xi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),mydata.yi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,1),mydata.vecImg(rngx,rngy,rngz,3),mydata.vecImg(rngx,rngy,rngz,2),2,'color','k'); %,0.08,'nointerp');
    hold off;
end
h(2)=subplot(2,2,2);
imagesc(squeeze(mydata.img(mydata.pos(1),:,:))'); d(2,:)=[2 3];
if ~isempty(mydata.vecImg)
    hold on;
    rngx=mydata.pos(1); rngy=1:5:size(mydata.xi,2); rngz=1:5:size(mydata.xi,3);
    quiver3(mydata.yi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),mydata.xi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,2),mydata.vecImg(rngx,rngy,rngz,3),mydata.vecImg(rngx,rngy,rngz,1),2,'color','k'); %,0.08,'nointerp');
    hold off;
end
h(3)=subplot(2,2,3);
imagesc(squeeze(mydata.img(:,:,mydata.pos(3)))'); d(3,:)=[1 2];
if ~isempty(mydata.vecImg)
    hold on;
    rngx=1:5:size(mydata.xi,1); rngy=1:5:size(mydata.xi,2); rngz=mydata.pos(3);
    quiver3(mydata.xi(rngx,rngy,rngz),mydata.yi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,1),mydata.vecImg(rngx,rngy,rngz,2),mydata.vecImg(rngx,rngy,rngz,3),2,'color','k'); %,0.08,'nointerp');
    hold off;
end

% some aesthetics
val = mydata.img(mydata.pos(1),mydata.pos(2),mydata.pos(3));
for i=1:3
    subplot(2,2,i);
    hold on; plot(mydata.pos(d(i,1)),mydata.pos(d(i,2)),'o','color',ones(1,3)*0.4,'linewidth',3,'markersize',12); hold off;
    axis xy; axis equal; axis tight; axis off; caxis(mydata.clim);
    title([num2str(mydata.pos(d(i,:))) ': ' num2str(val,'%.2f')])
end

if ~isempty(mydata.mri2mni)
    h(4) = subplot(2,2,4); axis off; caxis(mydata.clim);
    mniCoord = round(mydata.mri2mni*[mydata.pos 1]');
    coordInfo = {['Voxel: ' num2str(mydata.pos(1)) ',' num2str(mydata.pos(2)) ',' num2str(mydata.pos(3))],...
        ['MNI: ' num2str(mniCoord(1)) ',' num2str(mniCoord(2)) ',' num2str(mniCoord(3))]};
    title(h(4), coordInfo,'FontSize',16);
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