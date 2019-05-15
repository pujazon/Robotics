%Variables

global xmin;
global xmax;
global ymin;
global ymax;

%Left is ymin and top is xmin
xmin = 70;
xmax = 120;
ymin = 105;
ymax = 155;

%% Load

[VV,ss,dd] = dicomreadVolume(fullfile('./Skull_Tumor/Skull_Tumor/')); 
fv = stlread('Crani.stl');
NSlides = 112;
VV = squeeze(VV);
%Begins from 112 not from 1
%Work with two first in order to align
tumor_begin_slice = NSlides-40;
tumor_end_slice = NSlides-20;
%Color
colormap 'Bone';
tumor_slices = tumor_begin_slice:tumor_end_slice;

%% Keep only Tumor 

for k=1:size(VV,3)
    for i=1:size(VV,1)
        for j=1:size(VV,2)
            if(isTumor(i,j) == 0)
                VV(i,j,k)=0;
            end
        end
    end
end

%% Move Tumor display

xn = (-135:120)';
yn = (-150:105)';
zn = (1:112)';

contourslice(xn,yn,zn,VV,[],[],tumor_slices,2);
hold on

%% Translate Skull
T= transl(0,0,25);

for i=1:size(fv.vertices,1)
    aux = [fv.vertices(i,:) 1]';
    P = T*aux;
    fv.vertices(i,:) = P(1:3);
end

printSkull;

view(3);
axis tight

%%

function ret=isTumor(i,j)

    global xmin;
    global xmax;
    global ymin;
    global ymax;

    ret = 0;

    if((i>xmin)&&(i<xmax)&&(j>ymin)&&(j<ymax))
        ret = 1;
    end
   
end