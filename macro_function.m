%% macro spread simulation 1 may 23
function [area_bacteria_colony] = macro_function(dt,FluidPumpingRateFactor,imFlag,flag,FluidEffectConst,T,YeastDilutionFactor,GrowthrateMultiplier,AgnetcountMultiplier,YeastDivisionTime,BacteriaDivisionTime)
%% input variarables
count = 0;
% flag = 1; % 0==no yeast lawn; 1== with yeast and fluide; 2 == wiht yeast but no fluid; 3 = no pa
GriD = 0.01;
Max_nutrient = 10;% mol/kg

Yeast_inocu_conc = floor(YeastDilutionFactor*AgnetcountMultiplier*2.25e7); % number for cells 1 OD = 1.5e7 cells/ml $$ exp 500ul
Bact_inocu_conc = floor(AgnetcountMultiplier*2.4e6);  % number for cells 1OD = 8e8 cells/ml $$ exp 2ul

K2_nu = 1e-1; % nutrient weightage value = 0-1
K3_YeFluid = 1; % fluid weightage value = 0-1
kk_scaleConvert = 1; % 1 for cm
t_time_scale = 3600; % 1 for sec; 3600 for hrs

YeastConsuptionRate = t_time_scale*4.4e-12; % m/cell/sec
BacteriaConsuptionRate = t_time_scale*1.94e-12; % m/cell/sec
FluidPumpRate = FluidPumpingRateFactor*t_time_scale*0.3e-9; %cm3/sec
% FLuidEffect = 1e-6;% 1e-5 to 1e-7

D_nutr = kk_scaleConvert*kk_scaleConvert*t_time_scale*1.24e-5;%cm^2/sec
D_Fluid = kk_scaleConvert*kk_scaleConvert*t_time_scale*2.3e-5; %cm^2/sec

% dt = 1/t_time_scale; % i.e. 1e-4 hrs
dx = GriD;

lam = kk_scaleConvert*0.5e-4; % cm i.e. 0.5um
area_bacteria = pi*(lam*0.5)*(lam*1.5); % asuumed it to be a spherocylindrical
max_area = kk_scaleConvert*kk_scaleConvert*(10*GriD)*(10*GriD); %0.1 cm X 0.1 cm area of a box

Yeast_radi = kk_scaleConvert*4.5e-4; %cm i.e. 4.5um
area_yeast = pi*Yeast_radi*Yeast_radi; % asuumed it to be a sphere
NYeastSeed = 10000; % count of random yeast seed location on the grid

SeedPosition_Bacteria = round([0.5,0.5].*(1/GriD)); % bacteria will be in center always

if flag == 0 % i.e. if no yeast only bacteria
    seed_ya = zeros(1,NYeastSeed);
else % rand;; func of conc
    if NYeastSeed>Yeast_inocu_conc
        tem0 = rem(Yeast_inocu_conc,10);
        NYeastSeed = Yeast_inocu_conc - tem0;
    end

    seed_ya = diff([0,sort(randperm(Yeast_inocu_conc-1,NYeastSeed-1)),Yeast_inocu_conc]);
    % seed_ya = (Yeast_inocu_conc/NYeastSeed).*ones([1,NYeastSeed]);%round(0.05*(max_area/area_ye)); % i.e. initialize with 5% (hit and try) of maximum yeast capacity
end

if flag == 3 % i.e. if no bacteria only yeast lawn
    Seed_Bacteria = 0;
else
    Seed_Bacteria = Bact_inocu_conc;%round(0.9*(max_area/area_pa)); % i.e. initialize with 90% (hit and try) of maximum bacteria capacity
end

% Growth rate estimation
Yeast_Min_divisionTime = GrowthrateMultiplier*YeastDivisionTime(1)/t_time_scale;%sec 4500
Yeast_Max_divisionTime = GrowthrateMultiplier*YeastDivisionTime(2)/t_time_scale; %sec 4500
Bacteria_Min_divisionTime = GrowthrateMultiplier*BacteriaDivisionTime(1)/t_time_scale; %sec 2100
Bacteria_Max_divisionTime = GrowthrateMultiplier*BacteriaDivisionTime(2)/t_time_scale; %sec 2100

GrowthRateYeast = mean([Yeast_Max_divisionTime,Yeast_Min_divisionTime]);
GrowthRateYeast_diff = 0.5*(Yeast_Max_divisionTime-Yeast_Min_divisionTime);

GrowthRateBacteria = mean([Bacteria_Max_divisionTime,Bacteria_Min_divisionTime]);
GrowthRateBacteria_diff =0.5*(Bacteria_Max_divisionTime-Bacteria_Min_divisionTime);

%% initilization
xx  = GriD:GriD:1;
ImageBacteria = zeros(length(xx)); % name of varibale holding bacteria count on grid
ImageYeast = zeros([length(xx)]); % name of varibale holding yeast count on grid
ImageYeast_Fluid = zeros([length(xx)]); % name of varibale holding fluid on grid
Im_Nutrient = Max_nutrient.* ones(length(xx));% Initlizing nutrient on the grid
MAx_Area_grid = max_area.*ones([length(xx)]); % maximum area of each box on grid

Coordinate_rand_Loc_yeast = zeros([NYeastSeed,2]);

Im_boundary_condition = ones(length(xx)); %Dirichlet boundary condition
Im_boundary_condition(1,:) = 0;
Im_boundary_condition(end,:) = 0;
Im_boundary_condition(:,1) = 0;
Im_boundary_condition(:,end) = 0;

DiffusionEq_Sol_temp = zeros(length(xx)); % temporary derivatives variable

ImageBacteria(SeedPosition_Bacteria(1),SeedPosition_Bacteria(2)) = Seed_Bacteria; % Seeding bacteria

CarryCapacity_Bacteria = max_area.*ones(size(ImageBacteria)); %maximum carrying capacity of bacteria in terms of area
CarryCapacity_Yeast = max_area.*ones(size(ImageYeast)); %maximum carrying capacity of yeast in terms of area


Growth_Bacteria = zeros(size(ImageBacteria));
untruncated = makedist('Normal',GrowthRateBacteria,1);
if (GrowthRateBacteria-GrowthRateBacteria_diff) < (GrowthRateBacteria+GrowthRateBacteria_diff)
    truncated = truncate(untruncated,GrowthRateBacteria-GrowthRateBacteria_diff,GrowthRateBacteria+GrowthRateBacteria_diff);
    for i = 1:1/(GriD*GriD)
        r = random(truncated,1,1);
        Growth_Bacteria(i) = 1/r;
    end
else
    Growth_Bacteria =  (1/GrowthRateBacteria).*ones(size(ImageBacteria));
end

Growth_Yeast = zeros(size(ImageYeast));
untruncated = makedist('Normal',GrowthRateYeast,1);
if (GrowthRateYeast-GrowthRateYeast_diff) < (GrowthRateYeast+GrowthRateYeast_diff)
    truncated = truncate(untruncated,GrowthRateYeast-GrowthRateYeast_diff,GrowthRateYeast+GrowthRateYeast_diff);
    for i = 1:1/(GriD*GriD)
        r = random(truncated,1,1);
        Growth_Yeast(i) = 1/r;
    end
else
    Growth_Yeast =  (1/GrowthRateYeast).*ones(size(ImageYeast));
end

Overflow_temp_Bacteria = zeros(size(ImageBacteria));
Overflow_temp = zeros(size(ImageYeast));
% seeding Yeast
for i = 1:NYeastSeed
    Coordinate_rand_Loc_yeast(i,1) = ceil(rand(1)/GriD);
    Coordinate_rand_Loc_yeast(i,2) = ceil(rand(1)/GriD);
    ImageYeast(Coordinate_rand_Loc_yeast(i,1),Coordinate_rand_Loc_yeast(i,2)) = round(seed_ya(i)); % dont need rand seed ya
end
ImageYeast = ImageYeast.*Im_boundary_condition;

% neighbours list
% neighbours list
Val98 = length(xx);
TEmp_Neighbour = zeros([(Val98*Val98),8]); % list of all the possible 8 neighbours of a box on grid
for ii = 1:(length(xx)*length(xx))
    INdx = [ii-1,ii+1,ii-(Val98-1),ii+(Val98-1), ii-Val98,ii+Val98, ii-(Val98+1), ii+(Val98+1)];
    INdx(INdx<0 | INdx>(Val98*Val98)) = 0;
    TEmp_Neighbour(ii,:) = INdx; %boundar condition missing??
end

% neighbours counting yeast
Ya_layer = (ImageYeast>0); % converting yeast occupied region into binary image
ImageYeast_temp2 = eps + Ya_layer; % keep track of which box got filled first and which one was last

% neighbours counting bacteria
Bacteria_layer = (ImageBacteria>0); % converting bacteria occupied region into binary image
ImageBacteria_temp2 = eps + Bacteria_layer; % keep track of which box got filled first and which one was last

% Diffusion equation intialisation
dy = dx;
dd = dx*sqrt(2);
% 2D convolution masks - finite differences.
hN = [0 1 0; 0 -1 0; 0 0 0];
hS = [0 0 0; 0 -1 0; 0 1 0];
hE = [0 0 0; 0 -1 1; 0 0 0];
hW = [0 0 0; 1 -1 0; 0 0 0];
hNE = [0 0 1; 0 -1 0; 0 0 0];
hSE = [0 0 0; 0 -1 0; 0 0 1];
hSW = [0 0 0; 0 -1 0; 1 0 0];
hNW = [1 0 0; 0 -1 0; 0 0 0];

%% main {}
ik = 0;
YeastBact = Ya_layer + 5.*Bacteria_layer + ~Im_boundary_condition; % varible to keep track of free space
Time_tracker = dt;
while find(YeastBact ==0)>3
    ik = ik +1;

    % nutrient yeast and bacteria
    nablaN = imfilter(Im_Nutrient,hN,'conv','replicate');
    nablaS = imfilter(Im_Nutrient,hS,'conv','replicate');
    nablaW = imfilter(Im_Nutrient,hW,'conv','replicate');
    nablaE = imfilter(Im_Nutrient,hE,'conv','replicate');
    nablaNE = imfilter(Im_Nutrient,hNE,'conv','replicate');
    nablaSE = imfilter(Im_Nutrient,hSE,'conv','replicate');
    nablaSW = imfilter(Im_Nutrient,hSW,'conv','replicate');
    nablaNW = imfilter(Im_Nutrient,hNW,'conv','replicate');

    DiffusionEq_Sol_temp = D_nutr*(...
        (1/(dy^2))*nablaN + (1/(dy^2))*nablaS + ...
        (1/(dx^2))*nablaW + (1/(dx^2))*nablaE + ...
        (1/(dd^2))*nablaNE + (1/(dd^2))*nablaSE + ...
        (1/(dd^2))*nablaSW + (1/(dd^2))*nablaNW);

    % tt = imfilter(Im_Nutrient,kernelLApa,'replicate');
    if flag == 0 %no yeast lawn
        DiffusionEq_Sol_temp = dt.*(DiffusionEq_Sol_temp ...
            -(BacteriaConsuptionRate.*Bacteria_layer.*((ImageBacteria))));
    elseif flag == 3 %no bacteria
        DiffusionEq_Sol_temp = dt.*(DiffusionEq_Sol_temp...
            -(YeastConsuptionRate.*Ya_layer.*((ImageYeast))));
    else
        DiffusionEq_Sol_temp = dt.*(DiffusionEq_Sol_temp...
            -(YeastConsuptionRate.*Ya_layer.*((ImageYeast))) ...
            -(BacteriaConsuptionRate.*Bacteria_layer.*((ImageBacteria))));
    end
    Im_Nutrient = Im_Nutrient + DiffusionEq_Sol_temp;
    Im_Nutrient(Im_Nutrient<0) = 0; % in case, maximum nutrient is too low
    Im_Nutrient_probability = Im_Nutrient./Max_nutrient;

    %% fluid production and diffusion by yeast
    Temp_FluidProduc = (Ya_layer.*imcomplement(Bacteria_layer));

    nablaN = imfilter(ImageYeast_Fluid,hN,'conv','replicate');
    nablaS = imfilter(ImageYeast_Fluid,hS,'conv','replicate');
    nablaW = imfilter(ImageYeast_Fluid,hW,'conv','replicate');
    nablaE = imfilter(ImageYeast_Fluid,hE,'conv','replicate');
    nablaNE = imfilter(ImageYeast_Fluid,hNE,'conv','replicate');
    nablaSE = imfilter(ImageYeast_Fluid,hSE,'conv','replicate');
    nablaSW = imfilter(ImageYeast_Fluid,hSW,'conv','replicate');
    nablaNW = imfilter(ImageYeast_Fluid,hNW,'conv','replicate');

    DiffusionEq_Sol_temp = D_Fluid*(...
        (1/(dy^2))*nablaN + (1/(dy^2))*nablaS + ...
        (1/(dx^2))*nablaW + (1/(dx^2))*nablaE + ...
        (1/(dd^2))*nablaNE + (1/(dd^2))*nablaSE + ...
        (1/(dd^2))*nablaSW + (1/(dd^2))*nablaNW);


    DiffusionEq_Sol_temp = dt.*(DiffusionEq_Sol_temp...
        +(FluidPumpRate.*Temp_FluidProduc ...
        .*((ImageYeast.*Temp_FluidProduc))));

    ImageYeast_Fluid = ImageYeast_Fluid + DiffusionEq_Sol_temp;
    ImageYeast_Fluid_probability = (ImageYeast_Fluid./(eps+max(max(ImageYeast_Fluid))));

    if flag == 2 || flag == 0 %i.e. interaction without fluid || no yeast lawn
        ImageYeast_Fluid_probability = zeros(size(ImageBacteria));
        ImageYeast_Fluid = zeros(size(ImageBacteria));
    end

    % bacteria growth
    if flag == 1 % interaction with fluid film
        Im_combined_Bacteria = (K2_nu.*Im_Nutrient_probability + K3_YeFluid.* ImageYeast_Fluid_probability);
    else
        Im_combined_Bacteria = K2_nu.*Im_Nutrient_probability;
    end


    if flag == 2 || flag == 1 % growth of bacteria in presence of yeast
        CarryCapacity_Bacteria = (MAx_Area_grid - (ImageYeast.*area_yeast));
    end

    if flag ~=3 % i.e. in presence of bacteria
        ImageBacteria = ImageBacteria.*Im_boundary_condition;
        ImageBacteria = ImageBacteria + (ImageBacteria.*Growth_Bacteria.*Im_Nutrient_probability.*dt);
        ImageBacteria = ImageBacteria.*Im_boundary_condition;
        Fluid_importac = (ImageYeast_Fluid.* ImageYeast_Fluid)./(FluidEffectConst+(ImageYeast_Fluid.* ImageYeast_Fluid));
        Temp_Bacteria = floor((CarryCapacity_Bacteria.*(1-Fluid_importac))./area_bacteria); % maximum number of bacteria can be placed in a box
        Temp_Bacteria(Temp_Bacteria<0) = 0; %in case of overgrowth of yeast
        Overflow_temp_Bacteria = (floor(ImageBacteria) - Temp_Bacteria )>0; %8 is hit and try; below 8 we are getting infinite loop sometimes

        EXtra_Bacteria = floor(((ImageBacteria - Temp_Bacteria).*Overflow_temp_Bacteria)); % count of bacteria we need to relocate

        % Bacteria fountain distibution
        tt = (sum(Overflow_temp_Bacteria,'all'));
        Overflow_location_temp2 = find(Overflow_temp_Bacteria==1);
        [~,Overflow__temp] = sort(ImageBacteria_temp2(Overflow_location_temp2),'descend');
        Overflow_location_temp = Overflow_location_temp2(Overflow__temp);
        for i = 1:length(Overflow_location_temp)
            Overflow_loc = Overflow_location_temp(i);
            nEigh = TEmp_Neighbour(Overflow_loc,:);  % identify neighbours
            nEigh(nEigh==0) = []; % we are at the edge of the grid

            Differ_neigh =  ImageBacteria_temp2(nEigh); %Identify the iteraction at which box got filled
            nEigh(Differ_neigh>ImageBacteria_temp2(Overflow_loc)) = [];
            if isempty(nEigh) || length(nEigh)<3 % ideinfiy nighobur cells without backflow
                nEigh = TEmp_Neighbour(Overflow_loc,:);
                nEigh(nEigh==0) = [];
            end
            value = Im_combined_Bacteria(nEigh);
            RandomDist_weight = value./sum(value);
            Wight = round(EXtra_Bacteria(Overflow_loc).*RandomDist_weight);
            ImageBacteria(nEigh) = ImageBacteria(nEigh) + Wight;
            ImageBacteria(Overflow_loc) = ImageBacteria(Overflow_loc) - sum(round(Wight));
            Overflow_temp_Bacteria(Overflow_loc) = 0;
        end
        Bacteria_layer = (ImageBacteria>0);
        ImageBacteria_temp2 = ImageBacteria_temp2 + Bacteria_layer;
        ImageBacteria_temp2 = ImageBacteria_temp2.*Im_boundary_condition;
    end
    % yeast growth
    if flag ~=0  % i.e.with yeast lawn
        Im_combined_Ya = K2_nu.*Im_Nutrient_probability;
        if flag ~=3 % i.e.with bacteria
            CarryCapacity_Yeast = (MAx_Area_grid - (ImageBacteria.*area_bacteria));
        end
        ImageYeast = ImageYeast + (ImageYeast.*Growth_Yeast.*Im_Nutrient_probability.*dt);
        ImageYeast= ImageYeast.*Im_boundary_condition;
        Temp_ya = floor(CarryCapacity_Yeast/area_yeast);
        Temp_ya(Temp_ya<0) = 0;
        Overflow_temp = (floor(ImageYeast) - Temp_ya )>0;
        EXtra_ya = floor(((ImageYeast - Temp_ya).*Overflow_temp));
        %% Yeast fountain distibution
        Overflow_location_temp2 =find(Overflow_temp==1);
        [~,Overflow__temp] = sort(ImageYeast_temp2(Overflow_location_temp2),'descend');
        Overflow_location_temp = Overflow_location_temp2(Overflow__temp);
        for i = 1:length(Overflow_location_temp)
            Overflow_loc = Overflow_location_temp(i);
            nEigh = TEmp_Neighbour(Overflow_loc,:);
            nEigh(nEigh==0) = [];
            Differ_neigh =  ImageYeast_temp2(nEigh);
            nEigh(Differ_neigh>ImageYeast_temp2(Overflow_loc)) = [];
            if isempty(nEigh)|| length(nEigh)<3
                %                     Growth_YA(ji) = 0;
                nEigh = TEmp_Neighbour(Overflow_loc,:);
                nEigh(nEigh==0) = [];
            end
            value = Im_combined_Ya(nEigh);
            RandomDist_weight = value./sum(value);
            Wight = EXtra_ya(Overflow_loc).*RandomDist_weight;
            ImageYeast(nEigh) = ImageYeast(nEigh) + round(Wight);
            ImageYeast(Overflow_loc) = ImageYeast(Overflow_loc) - sum(round(Wight));
            Overflow_temp(Overflow_loc) = 0;
        end
        Ya_layer = (ImageYeast>0);
        ImageYeast_temp2 = ImageYeast_temp2 + Ya_layer;
        ImageYeast_temp2= ImageYeast_temp2.*Im_boundary_condition;
    end

    YeastBact = 2.*Ya_layer + 5.*Bacteria_layer + ~Im_boundary_condition;
    % Area of bacteria occupied region
    CC=bwconncomp(Bacteria_layer,8); %BW binary image ImSwarm_scale==1
    labeled=labelmatrix(CC);
    stats = regionprops(labeled,'Area');
    allAreas_Bacteria(ik) = sum([stats.Area]).*max_area;

    count = count +1;
    if imFlag == 1
        imageOut(max_area,area_yeast,area_bacteria,count,flag,ImageYeast,Ya_layer,ImageBacteria,Bacteria_layer,ImageBacteria_temp2,ImageYeast_temp2,ImageYeast_Fluid);
    end
    % break
    if Time_tracker>=T
        break;
    else
        Time_tracker = Time_tracker + dt;
    end
end
%% area of bacteria colony after 6 hours of simualtion
area_bacteria_colony = allAreas_Bacteria(ik);
end
