clear; close; clc; tic
%% input
NN = 10; % number of realizations
flag6 = zeros([NN,1]); % without yeast
flag6 = [flag6; 2*ones([NN,1])]; % with yeast without fluid pool
flag6 = [flag6; ones([NN,1])]; % with yeast with fluid pool
% flag6 = [flag6; 3*ones([NN,1])]; without bacteria
YeastDivisionTime = [4500;4500]; %sec %division time [min max]
BacteriaDivisionTime = [2100;2100]; %sec %division time [min max]

dt = 1/3600;
FluidPumpingRateFactor = 1;
FluidEffectConst = 0.9e-5;% a const
T = 6; % hours
GrowthrateMultiplier = 1;% cell growth rate ratio multiplier
AgnetcountMultiplier = 0.02;% cell count ratio multiplier
YeastDilutionFactor = 1;
imFlag = 0; % if want to extract each frames imFlage == 1

n = length(flag6);
area_pa = nan([size(flag6)]);
area_mean_temp = nan([size(YeastDivisionTime,2),length(YeastDilutionFactor),n ]);
area_mean_temp2 = nan([size(YeastDivisionTime,2),length(YeastDilutionFactor)]);

if imFlag == 1 % create direcotry for all the frames generated in each flag
    if not(isfolder('Frames'))
        mkdir('Frames');

        mkdir('Frames\YeastLayer');
        mkdir('Frames\YeastLayer\0');
        mkdir('Frames\YeastLayer\1');
        mkdir('Frames\YeastLayer\2');

        mkdir('Frames\CoExit');
        mkdir('Frames\CoExit\0');
        mkdir('Frames\CoExit\1');
        mkdir('Frames\CoExit\2');

        mkdir('Frames\FluidLayer');
        mkdir('Frames\FluidLayer\0');
        mkdir('Frames\FluidLayer\1');
        mkdir('Frames\FluidLayer\2');
    end
end

parfor ip = 1:n
    area_pa(ip) =  macro_function(dt,FluidPumpingRateFactor,imFlag,flag6(ip),FluidEffectConst,T,YeastDilutionFactor,GrowthrateMultiplier,AgnetcountMultiplier,YeastDivisionTime,BacteriaDivisionTime);
end

FlagType = unique(flag6);
if sum(FlagType==0)
    length3 = find(flag6 == 0);
    area_pa_0 = mean(area_pa(min(length3):max(length3)));
else
    area_pa_0 = 0;
end
if sum(FlagType==1)
    length3 = find(flag6 == 1);
    area_pa_1 = mean(area_pa(min(length3):max(length3)));
else
    area_pa_1 = 0;
end
if sum(FlagType==2)
    length3 = find(flag6 == 2);
    area_pa_2 = mean(area_pa(min(length3):max(length3)));
else
    area_pa_2 = 0;
end

X = categorical({'no Yeast','+Yeast -Fluid','+Yeast + Fluid'});
X = reordercats(X,{'no Yeast','+Yeast -Fluid','+Yeast + Fluid'});
Y = [area_pa_0, area_pa_2, area_pa_1];
bar(X,Y);
set(gca, 'FontSize',14)
ylabel('Simulated bacteria colony area (cm^{2})')

toc
