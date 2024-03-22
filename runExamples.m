clear;

roast
%pause
close all

roast('','','multipriors','on')
%pause
close all

roast([],{'Exx19',1,'C4',-1},'zeropadding',20)
%pause
close all

roast([],{'Exx19',1,'C4',-1},'zeropadding',20,'multipriors','on')
%pause
close all

roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4})
%pause
close all

roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'multipriors','on')
%pause
close all

roast('example/subject1.nii',[],'T2','example/subject1_T2.nii')
%pause
close all

roast('example/subject1.nii',[],'T2','example/subject1_T2.nii','multipriors','on')
%pause
close all

roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'resampling','on','zeropadding',20)
%pause
close all

roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'resampling','on','zeropadding',20,'multipriors','on')
%pause
close all

% roast('example/subject1.nii',{'G12',1,'J7',-1},'captype','biosemi')
% %pause
% close all

roast('example/subject1.nii',{'G12',0.25,'J7',-0.25,'Nk1',0.5,'Nk3',-0.5,'custom1',0.25,'custom3',-0.25},'captype','biosemi')
%pause
close all

roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'})
%pause
close all

roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'},'elecsize',{[8 2],[45 25 4],[5 8 2]})
%pause
close all

roast([],[],'electype','pad','elecori','ap')
%pause
close all

roast([],[],'electype','pad','elecori',[0.71 0.71 0])
%pause
close all

%roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',{'ap','lr','si'})
%pause
close all

roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',[0.71 0.71 0;-0.71 0.71 0;0 0.71 0.71])
%pause
close all % NOTE the direction vector for pads do not really work or intuitive

roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',[0.71 0.71 0;0 0.71 0.71])
%pause
close all

roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',{'ap',[],[0 0.71 0.71]})
%pause
close all

roast([],[],'meshoptions',struct('radbound',4,'maxvol',8))
%pause
close all

roast('example/subject1.nii',[],'resampling','on','zeropadding',10,'T2','example/subject1_T2.nii','electype','pad')
%pause
close all

roast([],{'Fp1',1,'FC4',1,'POz',-2},'conductivities',struct('csf',0.6,'electrode',0.1))
%pause
close all

roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'conductivities',struct('gel',[1 0.3 1],'electrode',[0.1 5.9e7 0.1]))
%pause
close all

roast('example/subject1.nii',{'Fp1',0.3,'F8',0.2,'POz',-0.4,'Nk1',0.5,'custom1',-0.6},...
        'electype',{'disc','ring','pad','ring','pad'},...
        'elecsize',{[],[7 9 3],[40 20 4],[],[]},...
        'elecori','ap','T2','example/subject1_T2.nii',...
        'meshoptions',struct('radbound',4,'maxvol',8),...
        'conductivities',struct('csf',0.6,'skin',1.0),...
        'resampling','on','zeropadding',30,...
        'simulationTag','awesomeSimulation')
%pause
close all

roast('nyhead')
%pause
close all

roast('nyhead',[],'multipriors','on')
%pause
close all
% doesn't skip all the way despite being the same
roast('nyhead',[],'resampling','on','zeropadding',25)
%pause
close all
% Error using roast
% The beauty of New York head is its 0.5 mm resolution. It's a bad practice to resample it
% into 1 mm. Use another head 'example/MNI152_T1_1mm.nii' for 1 mm model.
% 
% Error in runExamples (line 118)
roast('nyhead',[],'resampling','on','zeropadding',25,'multipriors','on')
%pause
close all

roast('nyhead',[],'electype','ring','elecsize',[7 10 3])
%pause
close all

roast('nyhead',{'Fp1',1,'FC4',1,'POz',-2},'electype','ring','elecsize',[7 10 3;6 8 3;4 6 2])
%pause
close all
