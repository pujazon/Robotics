%Variables

global xmin;
global xmax;
global ymin;
global ymax;

xmin = 100;
xmax = 150;
ymin = 70;
ymax = 120;

[VV,ss,dd] = dicomreadVolume(fullfile('./Skull_Tumor/Skull_Tumor/')); 
fv = stlread('Crani.stl');
VV = squeeze(VV);
scanner_position = ss.PatientPositions;
NSlides = 112;

%Begins from 112 not from 1
%Work with two first in order to align
tumor_begin_slice = NSlides-42;
tumor_end_slice = NSlides-24;

%Color
colormap 'Bone';
tumor_slices = tumor_begin_slice:tumor_end_slice;

xn = (-135:120)';
yn = (-150:105)';
zn = (1:112)';

contourslice(xn,yn,zn,VV,[],[],tumor_slices,2);
hold on

%Translate Skull
T= transl(0,0,25);
for i=1:size(fv.vertices,1)
    aux = [fv.vertices(i,:) 1]';
    P = T*aux;
    fv.vertices(i,:) = P(1:3);
end

printSkull;

view(3);
axis tight