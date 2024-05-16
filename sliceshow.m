function sliceshow(img,pos,color,clim,label,figName,vecImg,mri2mni,bbox)
% sliceshow displays a 3D volume allowing the user to click on different
% slices. If the variable "vecImg" is specified, sliceshow will also
% overlay a 3D vectorial field represented by arrows on the slices.
%
% If pos is given as input, then the function will start the display with
% the slices intersecting this position. Defaults to the middle of the
% volume.
% 
% Coordinates are editable. Enter the desired voxel or MNI coordinates to
% navigate to the desired position.
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
%
% bbox is the output of brainCrop(). Use this to display only the brain.
%
% sliceshow can display both the voxel and MNI coordinates.
%
% (c) Nov 02, 2017, Lucas C Parra
% (c) June 29, 2018, Yu (Andy) Huang
% (c) August 2019, Yu (Andy) Huang
% (c) 2024, Gavin Hsu and Andrew Birnbaum

if nargin >= 2 && ~isempty(bbox)
    minR = bbox(1,1);
    maxR = bbox(2,1);
    minA = bbox(1,2);
    maxA = bbox(2,2);
    minS = bbox(1,3);
    maxS = bbox(2,3);    
    img = img(minR:maxR, minA:maxA, minS:maxS);
else
    minR = 0;
    minA = 0;
    minS = 0;
end

if nargin<1 || isempty(img)
    error('At least give us a volume to display')
else
    [Nx,Ny,Nz] = size(img);
    mydata.img = img;
end

if nargin<2 || isempty(pos)
    mydata.pos = round(size(img)/2);
    mydata.pos_crop = mydata.pos + [minR,minA,minS];
else
     mydata.pos=pos-[minR,minA,minS];
     mydata.pos_crop = mydata.pos + [minR,minA,minS];
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
    if nargin >= 2 && ~isempty(bbox)
        vecImg = vecImg(minR:maxR, minA:maxA, minS:maxS,:);
    end
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
    mydata.mnipos = round(mydata.mri2mni*[mydata.pos_crop 1]');
end

whratio = 1.0187; %Width-to-height ratio
w = 8.5; %Width in inches
% link the callback function to new figure
fh = uifigure('WindowButtonDownFcn',@(src, event) myCallback(src, event, minR, minA, minS),'Name',figName,'NumberTitle','off','Units','inches','Position',[0,0,w,w/whratio],'AutoResizeChildren','off','SizeChangedFcn',@(src,event) figResize(src,event,minR,minA,minS));
movegui(fh,'center')
% store data as figure property
set(fh,'UserData',mydata);
set(fh,'HandleVisibility','on')

% initial display
showimages(minR,minA,minS,fh)
end

function myCallback(fh,~,minR, minA, minS)
% Selects new position depending on which subplot is clicked on.
% Updates display information in mydata.

mydata = get(fh,'UserData');
pos    = get(gca,'CurrentPoint');
pos    = round(pos(1,1:2));
switch gca
    case mydata.h(1); mydata.pos([1 3]) = pos;
    case mydata.h(2); mydata.pos([2 3]) = pos;
    case mydata.h(3); mydata.pos([1 2]) = pos;
end
mydata.pos_crop = mydata.pos + [minR,minA,minS];
mydata.mnipos = round(mydata.mri2mni*[mydata.pos_crop 1]');

% if new position is valid, update mydata as figure property and display
% otherwise do nothing
if ~sum( mydata.pos<1 | mydata.pos>size(mydata.img) )
    set(fh,'UserData',mydata);
    showimages(minR,minA,minS,fh)
end
end

function showimages(minR,minA,minS,fh)
% Displays the image based on current UserData.

% grab the stuff to display from figure
mydata  = get(fh,'UserData');

% Show the images and keep track of the axes that are generated
mytile = tiledlayout(fh,2,2,'Padding','tight','TileSpacing','none');

h(1) = nexttile(mytile,1); 
imagesc(h(1),squeeze(mydata.img(:,mydata.pos(2),:))'); d(1,:)=[1 3];
if ~isempty(mydata.vecImg)
    hold(h(1),'on');
    rngx=1:5:size(mydata.xi,1); rngy=mydata.pos(2); rngz=1:5:size(mydata.xi,3);
    quiver3(h(1),mydata.xi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),mydata.yi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,1),mydata.vecImg(rngx,rngy,rngz,3),mydata.vecImg(rngx,rngy,rngz,2),2,'color','k');
    hold(h(1),'off');
end
axtoolbar(h(1),{});

h(2) = nexttile(mytile,2); 
imagesc(h(2),squeeze(mydata.img(mydata.pos(1),:,:))'); d(2,:)=[2 3];
if ~isempty(mydata.vecImg)
    hold(h(2),'on');
    rngx=mydata.pos(1); rngy=1:5:size(mydata.xi,2); rngz=1:5:size(mydata.xi,3);
    quiver3(h(2),mydata.yi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),mydata.xi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,2),mydata.vecImg(rngx,rngy,rngz,3),mydata.vecImg(rngx,rngy,rngz,1),2,'color','k');
    hold(h(2),'off');
end
axtoolbar(h(2),{});

h(3) = nexttile(mytile,3); 
imagesc(h(3),squeeze(mydata.img(:,:,mydata.pos(3)))'); d(3,:)=[1 2];
if ~isempty(mydata.vecImg)
    hold(h(3),'on');
    rngx=1:5:size(mydata.xi,1); rngy=1:5:size(mydata.xi,2); rngz=mydata.pos(3);
    quiver3(h(3),mydata.xi(rngx,rngy,rngz),mydata.yi(rngx,rngy,rngz),mydata.zi(rngx,rngy,rngz),...
        mydata.vecImg(rngx,rngy,rngz,1),mydata.vecImg(rngx,rngy,rngz,2),mydata.vecImg(rngx,rngy,rngz,3),2,'color','k');
    hold(h(3),'off');
end
axtoolbar(h(3),{});

% Aesthetics
val = mydata.img(mydata.pos(1),mydata.pos(2),mydata.pos(3)); %EF value
for i=1:3
    han = nexttile(mytile,i);
    % Place position marker
    hold(han,'on'); plot(han,mydata.pos(d(i,1)),mydata.pos(d(i,2)),'o','color','m','linewidth',3,'markersize',12); hold(han,'off');
    axis(han,'xy'); axis(han,'equal'); axis(han,'tight'); axis(han,'off'); clim(han,mydata.clim);
    stackorder = get(h(i),'Children');
    if size(stackorder,1)>2
        set(h(i),'Children',[stackorder(1),stackorder(2),stackorder(3)]);
    end
    % Show crosshairs
    xline(han,mydata.pos(d(i,1))); yline(han,mydata.pos(d(i,2)));
end

% Set x- and y-lims equal in all views (for uniform scaling)
xlim(h,[0,max(size(mydata.img))])
ylim(h,[0,max(size(mydata.img))])

h(4) = nexttile(mytile,4); axis(h(4),'off'); box(h(4),'off'); clim(h(4),mydata.clim); axtoolbar(h(4),{});
mydata.pos_crop = mydata.pos + [minR,minA,minS];

% Arrange user input and coordinates in a UI Panel
ip = uipanel(fh,'BackgroundColor','white','AutoResizeChildren','off');
ip.Units = 'normalized';
ip.Position = [0.5,0,0.5,0.3];
ip.BorderType = 'none';
ip.Units = 'pixels';
panelsize = ip.Position;
ip.Units = 'normalized';
l = 0.28*panelsize(3); %left reference point
b = 0.59*panelsize(4); %bottom reference point
FSEF = 24; %Font size for displaying EF
FSc = 20; %Font size for coordinates
padding = 10;
%Display value:
uilabel(ip,'Position',[l-3.4*FSc,b-3*FSc,FSEF*12,FSEF+padding],'Text',strcat(mydata.label),'FontSize',FSEF,'FontWeight','bold','HorizontalAlignment','center');
uilabel(ip,'Position',[l-3.4*FSc,b-4.5*FSc,FSEF*12,FSEF+padding],'Text',num2str(val,'%.2f'),'FontSize',FSEF,'FontWeight','bold','HorizontalAlignment','center');
%Display voxel coordinates (callback below):
uilabel(ip,'Position',[l-3.4*FSc,b+2*FSc,3.1*FSc,FSc+padding],'Text','Voxel','FontSize',FSc,'FontWeight','bold');
efx = uieditfield(ip,'Position',[l+FSc,b+2*FSc,2.5*FSc,FSc+padding],'Tag','x','FontSize',FSc,'ValueChangedFcn',@(src,event) updateX(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
efx.Value = num2str(mydata.pos_crop(1));
uilabel(ip,'Position',[l,b+2*FSc,1.25*FSc,FSc+padding],'Text','X','FontSize',FSc);
efy = uieditfield(ip,'Position',[l+5.25*FSc,b+2*FSc,2.5*FSc,FSc+padding],'Tag','y','FontSize',FSc,'ValueChangedFcn',@(src,event) updateX(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
efy.Value = num2str(mydata.pos_crop(2));
uilabel(ip,'Position',[l+4.25*FSc,b+2*FSc,1.25*FSc,FSc+padding],'Text','Y','FontSize',FSc);
efz = uieditfield(ip,'Position',[l+9.5*FSc,b+2*FSc,2.5*FSc,FSc+padding],'Tag','z','FontSize',FSc,'ValueChangedFcn',@(src,event) updateX(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
efz.Value = num2str(mydata.pos_crop(3));
uilabel(ip,'Position',[l+8.5*FSc,b+2*FSc,1.25*FSc,FSc+padding],'Text','Z','FontSize',FSc);

%Display MNI coordinates (callback below):
if ~isempty(mydata.mri2mni)
    uilabel(ip,'Position',[l-3.4*FSc,b,3.1*FSc,FSc+padding],'Text','MNI','FontSize',FSc,'FontWeight','bold');
    efxm = uieditfield(ip,'Position',[l+FSc,b,2.5*FSc,FSc+padding],'Tag','x','FontSize',FSc,'ValueChangedFcn',@(src,event) updateXMNI(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
    efxm.Value = num2str(mydata.mnipos(1));
    uilabel(ip,'Position',[l,b,1.25*FSc,FSc+padding],'Text','X','FontSize',FSc);
    efym = uieditfield(ip,'Position',[l+5.25*FSc,b,2.5*FSc,FSc+padding],'Tag','y','FontSize',FSc,'ValueChangedFcn',@(src,event) updateXMNI(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
    efym.Value = num2str(mydata.mnipos(2));
    uilabel(ip,'Position',[l+4.25*FSc,b,1.25*FSc,FSc+padding],'Text','Y','FontSize',FSc);
    efzm = uieditfield(ip,'Position',[l+9.5*FSc,b,2.5*FSc,FSc+padding],'Tag','z','FontSize',FSc,'ValueChangedFcn',@(src,event) updateXMNI(src,event,minR,minA,minS,fh),'HorizontalAlignment','right','FontWeight','bold');
    efzm.Value = num2str(mydata.mnipos(3));
    uilabel(ip,'Position',[l+8.5*FSc,b,1.25*FSc,FSc+padding],'Text','Z','FontSize',FSc);
end

% Display colorbar
h(5) = colorbar(nexttile(mytile,2),'southoutside');
set(h(5),'YAxisLocation','bottom','FontSize',16);
h(5).Label.String = mydata.label;

colormap(fh,mydata.color);
set(fh,'color','w')

mydata.h=h;
set(fh,'UserData',mydata);
end

function updateX(src,event,minR, minA, minS,fh)
% Update figure to display location specified by voxel coordinate input
mydata = get(fh,'UserData');
if ~isempty(str2double(event.Value))
    switch src.Tag
        case 'x'
            mydata.pos(1) = round(str2double(event.Value))-minR;
            mydata.pos_crop(1) = str2double(event.Value);
        case 'y'
            mydata.pos(2) = round(str2double(event.Value))-minA;
            mydata.pos_crop(2) = str2double(event.Value);
        case 'z'
            mydata.pos(3) = round(str2double(event.Value))-minS;
            mydata.pos_crop(3) = str2double(event.Value);
    end
    mydata.mnipos = round(mydata.mri2mni*[mydata.pos_crop 1]');
end
if ~sum( mydata.pos<1 | mydata.pos>size(mydata.img) )
    set(fh,'UserData',mydata);
    showimages(minR,minA,minS,fh)
end
end

function updateXMNI(src,event,minR, minA, minS,fh)
% Update figure to display location specified by MNI coordinate input
mydata = get(fh,'UserData');
if ~isempty(str2double(event.Value))
    switch src.Tag
        case 'x'
            mydata.mnipos(1) = round(str2double(event.Value));
        case 'y'
            mydata.mnipos(2) = round(str2double(event.Value));
        case 'z'
            mydata.mnipos(3) = round(str2double(event.Value));
    end
    mydata.pos_crop = round(mydata.mri2mni\mydata.mnipos);
    mydata.pos(1) = mydata.pos_crop(1)-minR;
    mydata.pos(2) = mydata.pos_crop(2)-minA;
    mydata.pos(3) = mydata.pos_crop(3)-minS;
end
if ~sum( mydata.pos<1 | mydata.pos>size(mydata.img) )
    set(fh,'UserData',mydata);
    showimages(minR,minA,minS,fh)
end
end

function figResize(fh,~,minR,minA,minS)
    if ~isempty(fh.Children)
        showimages(minR,minA,minS,fh)
    end
end