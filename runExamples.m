clear;

try
    roast
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 1:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('','','multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 2:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/bikson.hdr')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 3:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/bikson.img')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 4:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Exx19',1,'C4',-1},'zeropadding',20)
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 5:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Exx19',1,'C4',-1},'zeropadding',20,'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 6:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4})
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 7:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 8:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',[],'T2','example/subject1_T2.nii')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 9:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',[],'T2','example/subject1_T2.nii','multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 10:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'resampling','on','zeropadding',20)
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 11:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4},'resampling','on','zeropadding',20,'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 12:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

% roast('example/subject1.nii',{'G12',1,'J7',-1},'captype','biosemi')
% %pause
% close all

try
    roast('example/subject1.nii',{'G12',0.25,'J7',-0.25,'Nk1',0.5,'Nk3',-0.5,'custom1',0.25,'custom3',-0.25},'captype','biosemi')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 13:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'})
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 14:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'},'elecsize',{[8 2],[45 25 4],[5 8 2]})
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 15:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],[],'electype','pad','elecori','ap')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 16:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],[],'electype','pad','elecori',[0.71 0.71 0])
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 17:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

%roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',{'ap','lr','si'})
%pause
close all

try
    roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',[0.71 0.71 0;-0.71 0.71 0;0 0.71 0.71])
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 18:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all % NOTE the direction vector for pads do not really work or intuitive

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',[0.71 0.71 0;0 0.71 0.71])
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 19:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',{'ap',[],[0 0.71 0.71]})
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 20:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],[],'meshoptions',struct('radbound',4,'maxvol',8))
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 21:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',[],'resampling','on','zeropadding',10,'T2','example/subject1_T2.nii','electype','pad')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 22:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'conductivities',struct('csf',0.6,'electrode',0.1))
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 23:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'conductivities',struct('gel',[1 0.3 1],'electrode',[0.1 5.9e7 0.1]))
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 24:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('example/subject1.nii',{'Fp1',0.3,'F8',0.2,'POz',-0.4,'Nk1',0.5,'custom1',-0.6},...
        'electype',{'disc','ring','pad','ring','pad'},...
        'elecsize',{[],[7 9 3],[40 20 4],[],[]},...
        'elecori','ap','T2','example/subject1_T2.nii',...
        'meshoptions',struct('radbound',4,'maxvol',8),...
        'conductivities',struct('csf',0.6,'skin',1.0),...
        'resampling','on','zeropadding',30,...
        'simulationTag','awesomeSimulation')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 25:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 26:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead',[],'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 27:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead',[],'resampling','on','zeropadding',25)
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 28:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead',[],'resampling','on','zeropadding',25,'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 29:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead',[],'zeropadding',5,'electype','ring','elecsize',[7 10 3],'multiaxial','on')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 30:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast('nyhead',{'Fp1',1,'FC4',1,'POz',-2},'electype','ring','elecsize',[7 10 3;6 8 3;4 6 2])
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 31:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
%pause
close all

try
    roast([],'leadfield','multiaxial','on','zeropadding',10,'electype','ring','simulationtag','LFwithMA')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 32:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
% pause
close all

try
    roast('nyhead','leadfield','multiaxial','on','zeropadding',10,'electype','ring','simulationtag','nyLFwithMA')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 33:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
% pause
close all

try
    roast_target([],'LFwithMA',[],'targetingtag','basic')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 34:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
% pause
close all

try
    roast_target([],'LFwithMA',[52 184 72;25 80 72;139 171 72],'coordType','voxel','optType','wls-l1','k',0.002,'targetradius',4,'targetingTag','3targets_radialIn')
catch ME
    fid=fopen('errLog.txt','a');
    fprintf(fid,'error at running cmd 35:\n');
    fprintf(fid,'%s\n',ME.message);
    fclose(fid);
end
% pause
close all

% add test on roast manual gui and align
% add test on reviewRes
