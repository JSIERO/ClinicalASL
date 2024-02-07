function SaveFIGUREtoPNG(data, mask, datarange, outputloc, suffix, label, COLORMAP, blackfont)
% ClinicalASL toolbox 2023, JCWSiero

% CBF or CVR map for different COLORMAPS: Crameri et al, "Misuse of colour in science communication' Nat. Comm 2020
% options are for example 'Viridis' (blue to yellow, converging map for CBF), 'vik' (Blue to Red) (VIK) (diverging map for CVR)
% please refrian from JET (non-scientific use only!)
switch COLORMAP
    case 'viridis' % now default for CBF maps: blue to yellow
        COLORMAP_black=viridis;
    case 'vik'  % now default for CVR maps: blue to white to red
        load('vik.mat');
        COLORMAP_black=vik;
    case 'jet'
        COLORMAP_black=jet;
    case 'imola'
        load('imola.mat');
        COLORMAP_black=imola;
    case 'parula'
        COLORMAP_black=parula;
    case 'devon'
        load('devon.mat');
        COLORMAP_black=devon;
    case 'turku'
        load('turku.mat');
        COLORMAP_black=turku;
end

COLORMAP_black(1,:)=[0 0 0];

if strcmp(label,'CBF')
    label2 = 'ml/100g/min';
elseif strcmp(label,'CVR')
    label2 = '\DeltaCBF ml/100g/min';
elseif strcmp(label,'time')
    label2 = 'time (s)';
elseif strcmp(label,'time_delta')
    label2 = '\Deltatime (s)';
elseif strcmp(label,'a.u.')
    label2 = 'a.u.';
elseif strcmp(label,'%')
    label2 = '%';
end
if isempty(mask)
    mask=ones(size(data));
end
data_masked=data.*mask;
data_masked(isnan(data_masked))=0;
image_data=immontage(rot90(data_masked),[ceil(size(data,3)/4),4],'noimage');

image_data((image_data<datarange(1) & image_data~=0))=0;
image_data(image_data>=datarange(2))=datarange(2);
Igz=find(image_data==0);

image_scaled=zeros(size(image_data));
image_scaled=255*(image_data - datarange(1))/(datarange(2) - datarange(1)) + 1;
image_scaled(Igz)=0;
RGB = ind2rgb(uint8(image_scaled),COLORMAP_black);
pause(0.5)
f1 = figure; 
figure(f1)
imagesc(RGB,datarange);
axis equal;
axis off;
set(f1, 'color', [0,0,0]); % black background color
colormap(COLORMAP_black);
h = colorbar('SouthOutside');
t=get(h,'Limits');
set(h,'Ticks',linspace(t(1),t(2),5))
set(get(h,'label'),'string',label2);

output = strcat(outputloc,suffix);
warning off
if nargin >=8
    if strcmp(blackfont,'yes')
        export_fig([output '_blkfnt'],'-png','-a4','-m4','-transparent','nocrop')
    end
end
h.Color=[0.85 0.85 0.85];
export_fig(output ,'-png','-a1','-m4','-nocrop')
warning on
close(f1)
end
