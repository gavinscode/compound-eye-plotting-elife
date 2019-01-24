%Written by Gavin Taylor, 2017. MIT License

%this script creates an arracy of structures containing information for each bee

%it is currently preloaded with the info. for the bees analyzed in Taylor et al. 2019
%which can be used as an example for other analysis

%Note that paths to downloaded data must be added where indicated by 'Fill here'

%FacetSizeFile is a csv file per bee that can be downloaded from doi:10.5061/dryad.23rj4pm
%the stackFolder is the eye label files which can be download from Morphosource* given the beeID
%the headStackFolder is the head label files which can be downloaded from Morphosource* given the headID

%*https://www.morphosource.org/Detail/ProjectDetail/Show/project_id/646
%reference the head/beeID to the scan# listed under 'element' in each media group for each sample

%analysis parameters for bee CE exterior analysis
clear; close all; clc

%parameters by bee
beeStructT = [];
beeStructT.BeeID = 'AM_60185';
beeStructT.BeeSpecies = 'AM';
beeStructT.BeeSex = 'F';
beeStructT.StackFolder = 'Fill here'; %add directory of label stack (and below)
    beeStructT.StackLabels.Lens = 2;
    beeStructT.StackLabels.Cones = 3;
    beeStructT.StackLabels.Retina = 4;
    beeStructT.StackLabels.LaminaC = 1;
    beeStructT.StackLabels.LensOuterS = 7;
    beeStructT.StackLabels.LensInnerS = 8;
    beeStructT.StackLabels.RetinaOuterS = 9;
    beeStructT.StackLabels.inside = 5;
beeStructT.LeftEyeTrans = reshape([1.69048 -0.299166 1.19636 0 0.088485 2.05415 0.388641 0 -1.23002 -0.263387 1.6722 0 6665.2 939.341 -334.077 1], [4,4]);
beeStructT.RightEyeTrans = reshape([1.00349 0 -0.0687345 0 0 1.00584 0 0 0.0687345 0 1.00349 0 -6415.85 -2.24243 507.921 1], [4,4]);
beeStructT.AddEyeScale = 2.09; %from left
beeStructT.VoxelSize = 5;
beeStructT.ITW = 2.9;
beeStructT.LeftEyeOriginal = 1;
beeStructT.headTrans = reshape([-0.175251 0.958411 0.225229 0 -0.978722 -0.144802 -0.145375 0 -0.106717 -0.245913 0.963389 0 9208.49 1369.26 -401.57 1], [4,4]);
beeStructT.headID = 'AM_60292';
beeStructT.headStackFolder = 'Fill here'; %add directory of head stack (and below)
beeStructT.headVoxelSize = 10; %actual scale should be five, but head resize takes care of this as x2 already in eye scale
beeStructT.FacetSizeFile = 'Fill here'; %add path to csv file (and below)
fullBeeData(1) = beeStructT;

beeStructT = [];
beeStructT.BeeID = 'AM_60186';
beeStructT.BeeSpecies = 'AM';
beeStructT.BeeSex = 'F';
beeStructT.StackFolder = 'Fill here';
    beeStructT.StackLabels.Lens = 2;
    beeStructT.StackLabels.Cones = 3;
    beeStructT.StackLabels.Retina = 4;
    beeStructT.StackLabels.LaminaC = 1;
    beeStructT.StackLabels.LensOuterS = 7;
    beeStructT.StackLabels.LensInnerS = 8;
    beeStructT.StackLabels.RetinaOuterS = 9;
    beeStructT.StackLabels.inside = 5;
beeStructT.LeftEyeTrans = reshape([0.0895934 -0.427387 -2.05036 0 -0.193453 -2.04513 0.417845 0 -2.08546 0.171349 -0.126853 0 7737.81 6172.16 2998.07 1], [4,4]);
beeStructT.RightEyeTrans = reshape([1 0 0 0 0 1 0 0 0 0 1 0 -5915.4 90.7571 -0.00012207 1], [4,4]);
beeStructT.AddEyeScale = 2.09; %from left
beeStructT.VoxelSize = 5;
beeStructT.ITW = 2.95;
beeStructT.LeftEyeOriginal = 1;
beeStructT.headTrans = reshape([-0.163008 0.965514 0.202998 0 -0.980478 -0.135618 -0.142299 0 -0.109864 -0.222234 0.968784 0 8227.32 865.35 -111.008 1], [4,4]);
beeStructT.headID = 'AM_60292';
beeStructT.headStackFolder = 'Fill here';
beeStructT.headVoxelSize = 10; %actual scale should be five, but head resize takes care of this as x2 already in eye scale
beeStructT.FacetSizeFile = 'Fill here';
fullBeeData(end+1) = beeStructT;

    beeStructT = [];
    beeStructT.BeeID = 'BT_77970';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([0.999403 0.031736 -0.0136017 0 -0.0313796 0.999178 0.0256578 0 0.0144048 -0.0252157 0.999579 0 3174.32 -10.798 -13.7057 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([-0.375712 -0.303546 0.958847 0 0.169199 -1.02805 -0.259158 0 0.99141 0.0604198 0.407595 0 225.074 3162.79 230.416 1], [4,4]);
    beeStructT.AddEyeScale = 1.07; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 4;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0514356 0.994221 0.0941647 0 -0.998222 -0.0540325 0.0252292 0 0.0301712 -0.0926999 0.995232 0 4405.25 -79.0945 -284.446 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;

    beeStructT = [];
    beeStructT.BeeID = 'BT_77971';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([0.99921 0.0379257 -0.0114999 0 -0.0373625 0.998248 0.0458048 0 0.0132171 -0.0453393 0.998884 0 3084.6 36.6367 -36.7945 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([-1.01321 -0.0148108 0.30406 0 -0.0355281 -1.04373 -0.16923 0 0.302339 -0.172284 0.999094 0 1061.15 3022.62 348.414 1], [4,4]);
    beeStructT.AddEyeScale = 1.05795; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 4.02;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0514356 0.994221 0.0941647 0 -0.998222 -0.0540325 0.0252292 0 0.0301712 -0.0926999 0.995232 0 4405.25 -79.0945 -284.446 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;

    beeStructT = [];
    beeStructT.BeeID = 'BT_77966';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([0.995691 0.0720628 0.0583422 0 -0.0740582 0.99671 0.0327962 0 -0.0557876 -0.0369756 0.997757 0 3181.66 -3.72021 -54.7035 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([0.129988 -0.0393972 -0.856649 0 -0.0337956 0.865522 -0.0449341 0 0.856887 0.0401127 0.128179 0 363.057 527.648 1197.2 1], [4,4]);
    beeStructT.AddEyeScale = 0.87; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 5.47;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0537743 0.993858 0.0967098 0 -0.998176 -0.0561625 0.0221459 0 0.0274411 -0.0953417 0.995056 0 4674.88 178.591 -269.511 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;

    beeStructT = [];
    beeStructT.BeeID = 'BT_77967';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([1.00285 0.0326124 -0.034362 0 -0.0318487 1.00322 0.0226436 0 0.0350723 -0.0215284 1.00313 0 3069.93 14.5735 7.86462 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([-0.783179 -0.00738203 0.320663 0 -0.0173652 -0.843871 -0.0618373 0 0.320276 -0.0638033 0.780761 0 1234.2 3025.73 280.022 1], [4,4]);
    beeStructT.AddEyeScale = 0.85; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 5.42;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0516125 0.998323 0.0260603 0 -0.99808 -0.0524624 0.0330842 0 0.0343962 -0.0243026 0.999111 0 4665.88 -87.6274 -134.487 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;
    
    beeStructT = [];
    beeStructT.BeeID = 'BT_77973';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([0.999282 0.0321818 -0.019837 0 -0.0315151 0.998958 0.0330465 0 0.0208807 -0.0323972 0.999251 0 2905.13 27.3152 -16.2327 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([-0.203248 -1.29794 0.231371 0 0.914247 -0.307425 -0.921486 0 0.949917 0.0181746 0.936384 0 160.585 3104.49 776.787 1], [4,4]);
    beeStructT.AddEyeScale = 1.33; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 1.97;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0508313 0.997896 0.0403382 0 -0.99801 -0.052268 0.0353451 0 0.0373789 -0.0384611 0.998563 0 4543.57 -90.561 -174.383 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;

    beeStructT = [];
    beeStructT.BeeID = 'BT_77974';
    beeStructT.BeeSpecies = 'BT';
    beeStructT.BeeSex = 'F';
    beeStructT.StackFolder = 'Fill here';
        beeStructT.StackLabels.Lens = 2;
        beeStructT.StackLabels.Cones = 3;
        beeStructT.StackLabels.Retina = 4;
        beeStructT.StackLabels.LaminaC = 1;
        beeStructT.StackLabels.LensOuterS = 7;
        beeStructT.StackLabels.LensInnerS = 8;
        beeStructT.StackLabels.RetinaOuterS = 9;
        beeStructT.StackLabels.inside = 5;
    beeStructT.LeftEyeTrans = reshape([0.998016 0.0565164 0.0276701 0 -0.0572728 0.997985 0.0273444 0 -0.0260688 -0.0288749 0.999241 0 2902.37 -10.916 -42.1242 1], [4,4]);
    beeStructT.RightEyeTrans = reshape([1.14563 -0.0193913 -0.472099 0 0.159503 1.18138 0.338539 0 0.444757 -0.37373 1.09463 0 -516.153 591.69 309.115 1], [4,4]);
    beeStructT.AddEyeScale = 1.24; %from right
    beeStructT.VoxelSize = 5;
    beeStructT.ITW = 2.97;
    beeStructT.LeftEyeOriginal = 0;
    beeStructT.headTrans = reshape([-0.0514356 0.994221 0.0941646 0 -0.998222 -0.0540325 0.0252292 0 0.0301712 -0.0926999 0.995232 0 4405.25 -79.0945 -284.446 1], [4,4]);
    beeStructT.headID = 'BT_60279';
    beeStructT.headStackFolder = 'Fill here';
    beeStructT.headVoxelSize = 5;
    beeStructT.FacetSizeFile = 'Fill here';
    fullBeeData(end+1) = beeStructT;
