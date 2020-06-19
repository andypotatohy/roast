% Maybe fieldtrip has to be updated?

%% Load the roast output and adapt it to brainstorm
%% load the headmodel
roast_headmodel = load('MNI152_T1_1mm_MNI152leadField.mat');
figure
for ind = 1 : 6%length(unique(roast_headmodel.elem(:,5)))
    % figure;
    hold on
    plotmesh(roast_headmodel.node(:,1:3), roast_headmodel.elem(roast_headmodel.elem(:,5)==ind,:));
    %   plotmesh(roast_headmodel.node(:,1:3), roast_headmodel.elem(roast_headmodel.elem(:,5)==ind,:),'facecolor','r');
end

% here you need to save the electrode position from the function  electrodePlacement.m line 158 
electrode_coord = load('lf_electrode_coord.mat');
hold on;
plotmesh(electrode_coord.electrode_coord,'ro')

%% Extract only the cortex 
cortexMeshElement = roast_headmodel.elem(find(roast_headmodel.elem(:,5)<=2),:);
figure;
plotmesh(roast_headmodel.node(:,1:3), cortexMeshElement,'x<50');

% [extractCortexSurface,elemid]  = volface(cortexMeshElement(:,1:4));
extractCortexSurface = volface(cortexMeshElement(:,1:4));
figure;
plotmesh(roast_headmodel.node(:,1:3), extractCortexSurface,'x<50');

surf.vertices = roast_headmodel.node(:,1:3);
surf.faces =extractCortexSurface;
 [surf,UsedV]=delete_unused_vertices(surf);
 surf.faces =     meshreorient( surf.vertices, surf.faces);
 figure;
 plotmesh(surf.vertices,surf.faces)

%% Load the leadfield
roast_leadfield = load('MNI152_T1_1mm_MNI152leadField_roastResult.mat');
roast_leadfield.A_all2 = permute(roast_leadfield.A_all,[3 1 2]);
roast_leadfield.A_cortex = roast_leadfield.A_all2(:,UsedV,:);
lfx = roast_leadfield.A_cortex(:,:,1); lfxVector = reshape(lfx,[],1);
lfy = roast_leadfield.A_cortex(:,:,2); lfyVector = reshape(lfy,[],1);
lfz = roast_leadfield.A_cortex(:,:,3); lfzVector = reshape(lfz,[],1);
lfAllVectors = [ lfxVector lfyVector lfzVector];
% matrice du gain
lfAll = reshape(lfAllVectors,size(roast_leadfield.A_cortex,1),[]);
% convert to brainstorm leadfield
bst_leadfield =[];
bst_leadfield.MEGMethod =  '';
bst_leadfield.EEGMethod =  'roast FEM';
bst_leadfield.ECOGMethod =  '';
bst_leadfield.SEEGMethod = '';
bst_leadfield.Gain = lfAll;%  [231�223983 double]
bst_leadfield.Comment = 'fromroast';
bst_leadfield.HeadModelType= 'surface';
bst_leadfield.GridLoc = surf.vertices;% [74661�3 double]
bst_leadfield.GridOrient = [];%: [74661�3 double]
bst_leadfield.GridAtlas = [];%
bst_leadfield.SurfaceFile = 'Subject01/tess_cortex_200607_2150.mat';
bst_leadfield.Param = [];
bst_leadfield.History = {'19-Jun-2015 00:08:01'  'import'  'Imported from Matlab variable: FEM'};

%% Save the cortex file to brainstorm
modelCortex.Comment = 'cortex roast';
modelCortex.Vertices = surf.vertices;
modelCortex.Faces = surf.faces;
modelCortex.VertConn = [];
modelCortex.VertNormals = [];
modelCortex.Curvature = [];
modelCortex.SulciMap = [];
% modelCortex.Atlas: [1�1 struct]
% modelCortex.iAtlas: 1
% modelCortex.tess2mri_interp: []
% modelCortex.Reg: []
% modelCortex.History: {2�3 cell}

%% Head model
% 1wm, 2gm, 3csf, 4skull, 5scalp, 6air
[no,el]=removeisolatednode(roast_headmodel.node(:,1:3),...
roast_headmodel.elem(roast_headmodel.elem(:,5)<=6,:));
newelem=meshreorient(no,el(:,1:4));
newelem = [newelem el(:,5)];
figure; plotmesh(no,el,'x>0')
% convert to brainstorm head model
FemMat = db_template('femmat');
FemMat.Comment = 'head model from roast';
FemMat.Vertices = no(:,1:3);
FemMat.Elements = newelem(:,1:4);
FemMat.Tissue = newelem(:,end);
FemMat.TissueLabels = [{'1-wm'} {'2-gm'} {'3-csf'} {'4-skull'} {'5-scalp'} {'6-air'} ];
FemMat.Tensors = [];
FemMat.History = [];


% %% Load the leadfield
% roast_leadfield = load('MNI152_T1_1mm_MNI152leadField_roastResult.mat');
% roast_leadfield.A_all2 = permute(roast_leadfield.A_all,[3 1 2]);
% lfx = roast_leadfield.A_all2(:,:,1); lfxVector = reshape(lfx,[],1);
% lfy = roast_leadfield.A_all2(:,:,2); lfyVector = reshape(lfy,[],1);
% lfz = roast_leadfield.A_all2(:,:,3); lfzVector = reshape(lfz,[],1);
% 
% lfAllVectors = [ lfxVector lfyVector lfzVector];
% lfAll = reshape(lfAllVectors,size(roast_leadfield.A_all,3),[]);

% A = [1:10; 11:20 ; 21:30]
% B= reshape(A,[],1)
% c =  reshape(B,3,[])


%% Load the channle file
bst_channel=[];
bst_channel.Comment =  ['roast channel (' num2str(length(electrode_coord.electrode_coord)) ')'];
bst_channel.MegRefCoef = [];
bst_channel.Projector = [];
bst_channel.TransfMeg = [];
bst_channel.TransfMegLabels = [];
bst_channel.TransfEeg = [];;
bst_channel.TransfEegLabels = {'manual correction'};
bst_channel.HeadPoints = [];
for iCh=1:length(electrode_coord.electrode_coord)-1
    bst_channel.Channel(iCh).Name = ['ch_' num2str(iCh)];
    bst_channel.Channel(iCh).Loc = electrode_coord.electrode_coord(iCh,:)';
    bst_channel.Channel(iCh).Orient = [];
    bst_channel.Channel(iCh).Comment = '';
    bst_channel.Channel(iCh).Type = 'EEG';
    bst_channel.Channel(iCh).Group = [];
end
bst_channel.IntraElectrodes = [];
bst_channel.History = channles.History;
bst_channel.SCS = channles.History;

% save('bst_channel','bst_channel')

