%% Analysis for for Duman & Azizi (2022)

%% Questions of Interest
% 1) Does elbow extension differ after nerve reinnervation? 
%    (extension rate, onset timing, & extension at touchdown)
% 2) Is hindlimb/ankle extension affected by nerve reinnervation? 
%    (extension rate, duration, & at touchdown)
% 3) Do activity of the Anconeus & Plataris change as a result of nerve
%    reinnervation? (timing of onset of plantaris prior to takeoff &
%    anconeus prior to landing/after takeoff)

% Read in Timing Data
Frames = readtable('Hindlimb_Proprioception_Timing.xlsx','Sheet',2,'Range','A:I');

% If inactive muscle times have already been saved then load them in
% if exist('t_muscle.mat') > 0
%     load t_muscle.mat
% end

for i = 524:size(Frames,1) % 1:size(Frames,1) to analyze all jumps
    % Calculating times of interest for jump
    times{i,1} = Find_Times(Frames,i);
    
    % Read Kinematic Data
    [Kin{i,1},L_hleg(i,1),trial] = Kinematics(Frames,i);

% QUESTION 1
    % Elbow Extension Rates
%     [Elbow_ang{i,1},dElbow_ang{i,1},avg_dElbExt(i,1),ElbExt_td(i,1),time_ElbExt_to(i,1),time_ElbExt_td(i,1)] = Calc_Elbow_Ext(Kin{i,1},times{i,1},i); % elbow (forelimb)

% QUESTION 2
    % Hindlimb & Ankle Extension Rates
%     [hLER{i,1},dhLER{i,1},avg_dhLER_to(i,1),max_dhLER_to(i,1),dhLER_to(i,1),t_HLflexes(i,1)] = Calc_hLER(Kin{i,1},times{i,1},L_hleg(i,1),i); % hindlimb

% QUESTION 3
    % Read & Analyze EMG Data (if exists)
    if i >= 301
        Folder = pwd; % get current folder's directory
        cd([Folder '\Igor Files']);
        Igor = readtable([trial '_emg.csv']);
        cd(Folder); % return to original directory
        
        E{i,1} = Trim_Igor(Igor,times{i,1});
%         if exist('t_muscle.mat') > 0 % use saved t_muscle values
%             [EMG_smooth{i,1},Plant_thresh{i,1},Anc_thresh{i,1},Plant_act_dur_to{i,1},Anc_ON_HLliftOff{i,1},Anc_act_dur_air{i,1},~] = Process_EMG(E{i,1},times{i,1},i,t_muscle{i,1});
%         else
            [EMG_smooth{i,1},Plant_thresh{i,1},Anc_thresh{i,1},Plant_act_dur_to{i,1},Anc_ON_HLliftOff{i,1},Anc_act_dur_air{i,1},t_muscle{i,1}] = Process_EMG(E{i,1},times{i,1},i);
%         end
        
    end

% Additional Variables
% Get the jump distance & treadmill belt movement
[belt_vel{i,1},avg_belt_vel(i,1),takeoff_vel(i,1),jump_dist(i,1)] = Calc_JumpParameters(Kin{i,1},times{i,1},i);

% Create column vector for toad number
toad(i,1) = str2num(trial(5:6));

% Display in Command Window that trial is done being analyzed
disp([Frames.trial{i} ' Jump: ' num2str(Frames.Jump(i)) ' Analyzed!']);
end

% Creating table for export (has same date, trial & jump # data as Frames)
Export_pt1 = Frames(:,1:3);

% Add variables of interst on to Export Table
% AAA_Export = [Export_pt1, table(toad,avg_dElbExt,ElbExt_td,time_ElbExt_to,time_ElbExt_td,avg_dhLER_to,max_dhLER_to,dhLER_to,t_HLflexes,avg_belt_vel,takeoff_vel,jump_dist,Plant_act_dur_to,Anc_ON_HLliftOff,Anc_act_dur_air)];
AAA_Export = [Export_pt1, table(toad,Plant_act_dur_to,Anc_ON_HLliftOff,Anc_act_dur_air)];
disp('Finished Analyzing ALL Trials!!!');

%% Saving Output Files
prompt = {"Would you like to save output table? Input Y or Yes to save, default setting is No.","Would you like to save muscle activation timing? Input Y or Yes to save, default setting is No."};
dlgtitle = 'Permission to Save Data';
dims = [1 55];
definput = {'No','No'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% Saving table
if strcmp(answer{1,1}(1),'Y')==1 || strcmp(answer{1,1}(1),'y')
    writetable(AAA_Export, 'Analyzed_Data.csv')
end

% Saving muscle inactive time values
if strcmp(answer{2,1}(1),'Y')==1 || strcmp(answer{2,1}(1),'y')
    save('t_muscle.mat','t_muscle')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 FUNCTIONS are BELOW in alphabetical order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calc_Elbow_Ext (function)
function [Elbow_ang,dElbow_ang,avg_dElbExt,ElbExt_td,time_ElbExt_to,time_ElbExt_td] = Calc_Elbow_Ext(K,t,row)
% function calculates the Elbow Angle (in degrees), its derivative, and 
% also pulls out points of interest like at touchdown.

%INPUTS
% K = kinematic points (points 5-8 are forelimb from metacarpophalangeal to
%     should joint, respectively)
% t = table containing times of interest where first row is relative to
%     hop and second row is relative to entire length of recording
% row = row number or trial number in Hindlimb_Proprioception_Timing
%       sequence

%OUTPUTS
% Elbow_ang = elbow angle in degrees over entire duration of trial
% dElbow_ang = derivative of elbow angle in degrees/second over entire duration of trial
% avg_ElbExt = average rate of elbow extension prior to landing
% ElbExt_td = elbow angle at touchdown (in degrees)
% time_ElbExt_to = time between hindlimb takeoff and beginning of elbow extension prior to landing
% time_ElbExt_td = time prior to touchdown when elbow first extends


% determining frames of interest
f = t; % setup table f (frames) to have same format as t (times)
if row < 301 % filmed at 100fps prior to September 25th 2020
    f{:,:} = round(f{:,:}*100 + 1);
    dt = 0.01; % time of one frame is 1/100th of a second
    numFrames = 1; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
else % filmed at 250fps from Septemeber 25th 2020 and after
    f{:,:} = round(f{:,:}*250 + 1);
    dt = 1/250; % time of 1 frame is 1/250th a second
    numFrames = 2; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
end

% creating vectors for each forelimb segment through duration of trial
for i = 1:size(K,1)
    v_fa(i,:) = K{i,4:6}-K{i,7:9}; % forearm
    u_fa(i,:) = v_fa(i,:)/sqrt(dot(v_fa(i,:),v_fa(i,:))); % unit direction of forearm away from elbow
    v_ua(i,:) = K{i,10:12}-K{i,7:9}; % upperarm
    u_ua(i,:) = v_ua(i,:)/sqrt(dot(v_ua(i,:),v_ua(i,:))); % unit direction of upperarm away from elbow
    Elbow_ang(i,1) = acosd(dot(u_ua(i,:),u_fa(i,:))); % elbow angle (degrees)
end
%Filter Elbow Angle with Savitzky-Golay Filter
if size(Elbow_ang,1) >=49
    Elbow_ang(:,1) = sgolayfilt(Elbow_ang(:,1),7,49);
elseif mod(size(Elbow_ang,1),2) == 1 % if remainder of the length divided by 2 is 1 (length of Elbow Angle is odd)
    Elbow_ang(:,1) = sgolayfilt(Elbow_ang(:,1),7,size(Elbow_ang,1));
else
    Elbow_ang(:,1) = sgolayfilt(Elbow_ang(:,1),7,(size(Elbow_ang,1)-1));
end

ElbExt_td = Elbow_ang(f{1,4},1);    % Elbow Extension at touchdown

dElbow_ang = zeros(length(Elbow_ang),1); % intial derivative of Elbow Angle (assumes stationary start)
for i = (numFrames+1):(length(Elbow_ang)-numFrames)
    dElbow_ang(i,1) = (Elbow_ang(i+numFrames,1)-Elbow_ang(i-numFrames,1))/(2*numFrames*dt); % rate of elbow flexion(-)/extension(+)
end
for i =  0:(numFrames-1)
    dElbow_ang(end-i,1) = dElbow_ang(end-numFrames,1); % set last few points to same as final calculated
end

% Calculating Average Rate of Elbow Extension prior to landing
dElbExt_sign = sign(dElbow_ang);
dElbExt_air = nan(size(dElbow_ang)); % define rate of elbow extension in air as NaNs to start (only care about extension just prior to landing)
n = size(dElbExt_air,1);
while dElbExt_sign(end,1) <= 0
    dElbExt_sign(end,:) = []; % delete instance at end when elbow is flexing (after touchdown)
    n = n - 1; % subtract off 1 to get next row of interest
end
while dElbExt_sign(end,1) > 0
    dElbExt_sign(end,:) = []; % delete instance where sign is positive (elbow is extending in air)
    dElbExt_air(n,1) = dElbow_ang(n,1); % saves rate of elbow extension to aerial phase vector
    n = n - 1; % subtract off 1 to get next row of interest
end
avg_dElbExt = mean(dElbExt_air,'omitnan');

% Calculating onset time of elbow extension relative to HL takeoff
Ext = dElbExt_air; % times when elbow is extending just prior to landing
while isnan(Ext(end,1))==1
    Ext(end,:) = []; % deletes last row while elbow is not extending
end
while isnan(Ext(end-1,1))==0
    Ext(end,:) = []; % deletes last row when 2nd to last row is still during elbow extension
end
if size(Elbow_ang,1) > f{1,3}
    time_ElbExt_to = (size(Ext,1)-f{1,3})*dt; % time between hindlimb takeoff and beginning of elbow extension prior to landing
else
    time_ElbExt_to = nan;
end
time_ElbExt_td = -(f{1,5}-size(Ext,1))*dt; % time prior to touchdown when elbow first extends
end

%% Calc_hLER (function)
function [hLER,dhLER,avg_dhLER_to,max_dhLER_to,dhLER_to,t_HLflexes] = Calc_hLER(K,t,L,row)
% function calculates the HindLimb Extension Ratio (hLER) of the HINDLIMB, its
% derivative (dhLER) and also pulls out the change in hLER during takeoff.

%INPUTS
% K = kinematic points (points 11 & 14 are hind-toe and COM estimate (posterior end of spine), respectively)
% t = table containing times of interest where first row is relative to
%     hop and second row is relative to entire length of recording
% L = hindlimb leg length

%OUTPUTS
% hLER = hindlimb extension ratio over entire duration of trial
% dhLER = derivative of forelimb extension ratio over duration of trial
% avg_dhLER_to = average rate of hindlimb extension during takeoff
% max_dhLER_to = maximum rate of hindlimb extension during takeoff
% dhLER_to = rate of hindlimb extension at takeoff
% t_HLflexes = time in ms when hindlimb extension ends and flexion begins

% determining frames of interest
f = t; % setup table f (frames) to have same format as t (times)
if row < 301 % filmed at 100fps prior to September 25th 2020
    f{:,:} = round(f{:,:}*100 + 1);
    dt = 1/100; % time in seconds between frames
    numFrames = 1; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
else % filmed at 250fps from Septemeber 25th 2020 and after
    f{:,:} = round(f{:,:}*250 + 1);
    dt = 1/250;
    numFrames = 2; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
end


% creating vector for hindlimb through duration of trial
v_hlimb = [K.HIx, K.HIy, K.HIz] - [K.TOx, K.TOy, K.TOz]; % whole hindlimb from toe (TO) to ischium (HI; hip estimate)
for i = 1:size(K,1)
    L_hlimb(i,1) = sqrt(dot(v_hlimb(i,:),v_hlimb(i,:))); % hindlimb length
    hLER(i,1) = L_hlimb(i,1)/L; % Limb Extension Ratio
end
% Filter hLER with Savitzky-Golay Filter
if size(hLER,1) >=49
    hLER(:,1) = sgolayfilt(hLER(:,1),7,49);
elseif mod(size(hLER,1),2) == 1 % if remainder of the length divided by 2 is 1 (length of LER is odd)
    hLER(:,1) = sgolayfilt(hLER(:,1),7,size(hLER,1));
else
    hLER(:,1) = sgolayfilt(hLER(:,1),7,(size(hLER,1)-1));
end

dhLER = zeros(length(hLER),1); % creates inital derivative of hLER vector, with all zeros (assumes stationary start)
for i = (numFrames+1):(length(hLER)-numFrames)
    dhLER(i,1) = (hLER(i+numFrames)-hLER(i-numFrames))/(2*numFrames*dt); % derivative of Limb Extesion Ratio
end
for i =  0:(numFrames-1)
    dhLER(end-i,1) = dhLER(end-numFrames,1); % assumes low acceleration in LER toward end of landing
end

if size(dhLER,1) < f{1,3}
    dhLER_to = nan;
    max_dhLER_to = nan;
    avg_dhLER_to = nan;
else
    dhLER_to = dhLER(f{1,3},1); % hindlimb extension rate at takeoff
    max_dhLER_to = max(dhLER(1:f{1,3},1)); % maximum rate of hindlimb extension prior to takeoff
    avg_dhLER_to = mean(dhLER(1:f{1,3},1)); % average rate of hindlimb extension prior to takeoff
end
if max_dhLER_to == 0 && isnan(dhLER_to) == 1 && isnan(avg_dhLER_to) == 1
    max_dhLER_to = nan; % reset to nan if all values are nan
end

% Calculating time at which Hindlimb Extension Stops during early part of
% hop (same as when HL flexion starts occuring)
HL_ext = dhLER;
if sum(HL_ext < 0)==0 % if Hindlimb is always extending
    t_HLflexes = nan;
else
    while HL_ext(end,1) >= 0
        HL_ext(end,:) = []; % deletes frames at end when hindlimb may be extending
    end
    while HL_ext(end,1) < 0
        HL_ext(end,:) = []; % deletes frames at end when hindlimb is flexing
    end
    t_HLflexes = length(HL_ext)*dt; % time in s when hindlimb extension ends and flexion begins
end

end

%% Calc_JumpParameters (function)
function [belt_vel,avg_belt_vel,takeoff_vel,jump_dist] = Calc_JumpParameters(K,t,row)
% Calculates the jump distance for a hop and the belt velocity through the
% entire time of interest

%INPUTS
% K = kinematics table
% t = times of interest where first row is in reference to jump and second
%     row is relative to entire 60s recording

%OUTPUTS
% belt_vel  = speed of belt over the duration of the jump trial of interest
% avg_belt_vel  = average speed of belt during entire jump trial (m/s)
% takeoff_vel = instantaneous velocity of the center of mass at hindlimb liftoff 
% jump_dist = distnace toad jumps based on movement of wrist from liftoff
%             to touchdown & adding in belt movement during that time.
%                 hindlimb liftoff


% determining frames of interest
f = t; % setup table f (frames) to have same format as t (times)
if row < 301 % filmed at 100fps prior to September 25th 2020
    f{:,:} = round(f{:,:}*100 + 1);
    dt = 1/100; % differential unit of time between frames
    numFrames = 1; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
else % filmed at 250fps from Septemeber 25th 2020 and after
    f{:,:} = round(f{:,:}*250 + 1);
    dt = 1/250; % differential unit of time between frames
    numFrames = 2; % number of frames to check over to calculate derivative (final - initial = N+numFrames - N-numFrames)
end

% Filter Kinematics Data
if size(K,1) >=49
    for i = 1:size(K,2)
    K_filt(:,i) = sgolayfilt(K{:,i},7,49); % Savitzky-Golay filter
    end
elseif mod(size(K,1),2) == 1 % if remainder of the length divided by 2 is 1 (length of data is odd)
    for i = 1:size(K,2)
    K_filt(:,i) = sgolayfilt(K{:,i},7,size(K,1));
    end
else
    for i = 1:size(K,2)
    K_filt(:,i) = sgolayfilt(K{:,i},7,(size(K,1)-1));
    end
end



% Treadmill Belt Velocity (meters per second)
belt_vel = zeros(size(K_filt,1),1); % define size & set all values to zero
for i = (numFrames+1):(size(K_filt,1)-numFrames)
    belt_vel(i,1) = (K_filt(i+numFrames,1)-K_filt(i-numFrames,1))/(2*numFrames*dt);
end
belt_vel(1,1) = belt_vel(2,1); % assume constant velocity
for i = 0:(numFrames-1)
    belt_vel(end-i,1) = belt_vel(end-numFrames,1); % assume constant velocity
end
avg_belt_vel = mean(belt_vel);

% Jump Distance (meters)
D_forelimb = K_filt(f{1,4},4:6) - K_filt(f{1,2},4:6); % displacement of wrist during jump
D_belt = K_filt(f{1,4},1:3) - K_filt(f{1,2},1:3); % displacement of treadmill belt during jump
jump_dist = sqrt(dot(D_forelimb,D_forelimb)) + sqrt(dot(D_belt,D_belt));

% COM Velocity (meters per second)
COM_disp = zeros(size(K_filt,1),3); % define size & set all values to zero
COM_vel = zeros(size(K_filt,1),1); % define size & set all values to zero
for i = 2:(size(K_filt,1)-1)
    COM_disp(i,:) = K_filt(i+1,25:27)-K_filt(i-1,25:27);
    COM_vel(i,1) = sqrt(dot(COM_disp(i,1),COM_disp(i,1)))/(dt*2);
end
COM_vel(1,1) = COM_vel(2,1)-(COM_vel(3,1)-COM_vel(2,1));
COM_vel(end,1) = COM_vel(end-1,1) + (COM_vel(end-1,1)-COM_vel(end-2,1));
if size(COM_vel,1) > f{1,3}
    takeoff_vel = COM_vel(f{1,3},1);
else
    takeoff_vel = nan;
end

end

%% Find_Times (function)
function t = Find_Times(F,row)
% This function finds the times (t) associated with the particular trial

% INPUTS
% F - the frames of interest for all jumps
% row - the row within frames for the jump trial of interest

% OUTPUTS
% t - table containing times of interest where first row is relative to
%     hop or trial and second row is relative to entire length of recording

if row < 301 % filmed prior to Septmeber 25th, 2020 (at 100fps)
    Start = table2array(F(row,4))/100; % times are in second since filmed at 100fps
    ForelimbUp = table2array(F(row,5))/100;
    HindlimbUp = table2array(F(row,6))/100;
    ForelimbDown = table2array(F(row,7))/100;
    BodyDown = table2array(F(row,8))/100;
else % filmed on or after September 25th, 2020 (at 250fps)
    Start = table2array(F(row,4))/250; % times are in second since filmed at 250fps
    ForelimbUp = table2array(F(row,5))/250;
    HindlimbUp = table2array(F(row,6))/250;
    ForelimbDown = table2array(F(row,7))/250;
    BodyDown = table2array(F(row,8))/250;
end
    

t = table(Start,ForelimbUp,HindlimbUp,ForelimbDown,BodyDown);
t(2,:) = t(1,:); % copies times to second row
t{1,:} = t{1,:}-(t{2,1});
end

%% Kinematics (function)
function [K,L_hleg,trial] = Kinematics(Frames,i)
% Reads in Kinematics Data

% INPUTS
% Frames = list of jump data
% i = row number for jump of interest

% OUTPUTS
% K - table of digitized kinematics points
% L_hleg - length of hindlimb
% trial - charater string of the trial name of interest

    % Initial folder
    Folder = pwd;

    % Reading in Kinematics Data
    trial = cell2mat(Frames{i,2});
    date = cell2mat(Frames{i,1});
    cd([Folder '\Digitized Files'])
    Digit = readtable([trial '_xyzpts.csv']);
    cd(Folder)
    
    if i < 301 % prior to September 25th, 2020
        % Naming Points
        Digit.Properties.VariableNames{1} = 'GNDx'; % ground or treadmill belt
        Digit.Properties.VariableNames{2} = 'GNDy';
        Digit.Properties.VariableNames{3} = 'GNDz';
        Digit.Properties.VariableNames{4} = 'WRx'; % wrist
        Digit.Properties.VariableNames{5} = 'WRy';
        Digit.Properties.VariableNames{6} = 'WRz';
        Digit.Properties.VariableNames{7} = 'ELx'; % elbow
        Digit.Properties.VariableNames{8} = 'ELy';
        Digit.Properties.VariableNames{9} = 'ELz';
        Digit.Properties.VariableNames{10} = 'SHx'; % shoulder
        Digit.Properties.VariableNames{11} = 'SHy';
        Digit.Properties.VariableNames{12} = 'SHz';
        Digit.Properties.VariableNames{13} = 'TOx'; % toe
        Digit.Properties.VariableNames{14} = 'TOy';
        Digit.Properties.VariableNames{15} = 'TOz';
        Digit.Properties.VariableNames{16} = 'MEx'; % metatarsals
        Digit.Properties.VariableNames{17} = 'MEy';
        Digit.Properties.VariableNames{18} = 'MEz';
        Digit.Properties.VariableNames{19} = 'ANx'; % ankle
        Digit.Properties.VariableNames{20} = 'ANy';
        Digit.Properties.VariableNames{21} = 'ANz';
        Digit.Properties.VariableNames{22} = 'KNx'; % knee
        Digit.Properties.VariableNames{23} = 'KNy';
        Digit.Properties.VariableNames{24} = 'KNz';
        Digit.Properties.VariableNames{25} = 'HIx'; % hip
        Digit.Properties.VariableNames{26} = 'HIy';
        Digit.Properties.VariableNames{27} = 'HIz';
    else % September 25th, 2020 and later
        cd(['C:\Users\AJD44\Desktop\Chapter 3\' date '\Toe Hip'])
        Digit_TH = readtable([trial '_TH_xyzpts.csv']);
        % Naming Points 
        Digit.Properties.VariableNames{1} = 'GNDx'; % ground or treadmill belt
        Digit.Properties.VariableNames{2} = 'GNDy';
        Digit.Properties.VariableNames{3} = 'GNDz';
        Digit.Properties.VariableNames{4} = 'WRx'; % wrist
        Digit.Properties.VariableNames{5} = 'WRy';
        Digit.Properties.VariableNames{6} = 'WRz';
        Digit.Properties.VariableNames{7} = 'ELx'; % elbow
        Digit.Properties.VariableNames{8} = 'ELy';
        Digit.Properties.VariableNames{9} = 'ELz';
        Digit.Properties.VariableNames{10} = 'SHx'; % shoulder
        Digit.Properties.VariableNames{11} = 'SHy';
        Digit.Properties.VariableNames{12} = 'SHz';
        Digit.Properties.VariableNames{13} = 'MEx'; % metatarsals
        Digit.Properties.VariableNames{14} = 'MEy';
        Digit.Properties.VariableNames{15} = 'MEz';
        Digit.Properties.VariableNames{16} = 'ANx'; % ankle
        Digit.Properties.VariableNames{17} = 'ANy';
        Digit.Properties.VariableNames{18} = 'ANz';
        Digit.Properties.VariableNames{19} = 'KNx'; % knee
        Digit.Properties.VariableNames{20} = 'KNy';
        Digit.Properties.VariableNames{21} = 'KNz';
        Digit_TH.Properties.VariableNames{1} = 'TOx'; % toe
        Digit_TH.Properties.VariableNames{2} = 'TOy';
        Digit_TH.Properties.VariableNames{3} = 'TOz';
        Digit_TH.Properties.VariableNames{4} = 'HIx'; % hip
        Digit_TH.Properties.VariableNames{5} = 'HIy';
        Digit_TH.Properties.VariableNames{6} = 'HIz';
        % Arranging to be same as previous trials
        Digit = [Digit(:,1:12),Digit_TH(:,1:3),Digit(:,13:21),Digit_TH(:,4:6)];
    end
    %Return to initial folder
    cd(Folder);

% Calculating time range of interest for specific trial
trial = char(Frames{i,2});
date = char(Frames{i,1});
jump = Frames{i,3};
% Removes trials from Frames that are not from same 60s video
while strcmp(trial,char(Frames{1,2}))~=1 || strcmp(date,char(Frames{1,1}))~=1
    Frames(1,:)=[];
end
while strcmp(trial,char(Frames{end,2}))~=1 || strcmp(date,char(Frames{end,1}))~=1
    Frames(end,:)=[];
end
% Determining absolute frame number in trimmed video (& kinematics files)
frames = NaN(size(Frames,1),2);

frames(1,1) = 1; % initial frame
L = Frames.BodyTouchdown(1)-Frames.Start(1);
frames(1,2) = frames(1,1)+L;

if size(Frames,1) > 1
    for j = 2:size(Frames,1)
        frames(j,1) = frames(j-1,2) + 1;
        L = Frames.BodyTouchdown(j) - Frames.Start(j);
        frames(j,2) = frames(j,1) + L;
    end
end
    
% Kinematics during jump of interest
K = Digit(frames(jump,1):frames(jump,2),:);

% determine hindlimb length
for j = 1:size(Digit,1)
    % leg segment vectors
    v_toe(j,:) = Digit{j,13:15} - Digit{j,16:18};
    v_foot(j,:) = Digit{j,16:18} - Digit{j,19:21};
    v_shank(j,:) = Digit{j,19:21} - Digit{j,22:24};
    v_thigh(j,:) = Digit{j,22:24} - Digit{j,25:27};
    % segment lengths
    L_toe(j,1) = sqrt(dot(v_toe(j,:),v_toe(j,:)));
    L_foot(j,1) = sqrt(dot(v_foot(j,:),v_foot(j,:)));
    L_shank(j,1) = sqrt(dot(v_shank(j,:),v_shank(j,:)));
    L_thigh(j,1) = sqrt(dot(v_thigh(j,:),v_thigh(j,:)));
end
L_hleg = mean(L_toe) + mean(L_foot) + mean(L_shank) + mean(L_thigh);
end

%% Process_EMG (function)
function [EMG_smooth,Plant_thresh,Anc_thresh,Plant_act_dur_to,Anc_ON_HLliftOff,Anc_act_dur_air,t_muscles] = Process_EMG(E,times,row,varargin)
% Function calculates the rectified and smoothed EMG traces for both
% muscles and also provides the threshold for activation for each.


% INPUTS
% E - table containing raw EMG traces for both muscles and starts 0.5
%     seconds prior to start of motion
% times - table containing times of interest where first row is relative to
%         hop and second row is relative to entire length of recording
% row - row or trial of interest
% varargin - t_muscles (optional variable) gives the frame numbers used 
%            during trial when muscles were inactive and active


% OUTPUTS
% EMG_smoothed - table of rectified and smoothed EMG traces
% Plant_thresh - activation threshold value for plantaris
% Anc_thresh - activation threshold value for anconeus
% Plant_act_dur_to - duration of time in seconds plantaris is active during
%                    takeoff
% Anc_ON_HLliftOff - time in seconds anconeus activates after takeoff
% Anc_act_dur_air - duration of time in seconds anconeus is active during
%                   aerial phase of hop
% t_muscles - gives the frame numbers used during trial when muscles were
%              inactive and active


% FORMER OUTPUTS
% Plant_ON_takeoff - time in seconds plantaris activates prior to takeoff
% Plant_ON_start - time in seconds plantaris activates relative to start of
%                  motion (trial)

%% Rectify & Smooth EMG signals
EMG_rect = abs(table2array(E)); % rectify traces
for i = 1:size(EMG_rect,2)
    EMG_smooth(:,i) = movmean(EMG_rect(:,i),9); % smoothing average
end

%% Find ranges where muscles are inactive and active
if isempty(varargin) == 1 % don't have activation timing so need to collect it
    R = num2str(row); % convert row or trial number to string
    figName = [R ': Select ranges of Plantaris, inactive then active'];
    figure('Name',figName)
    plot(E.plantaris,'k-')
    hold on;
    plot(EMG_smooth(:,1),'r-')
    xline(5001,'k--','Label','Start Moving');
    xline((times{1,2}*10000 + 5000),'Label','Forelimb Liftoff');
    xline((times{1,3}*10000 + 5000),'b-','Label','Hindlimb Liftoff');
    xline((times{1,4}*10000 + 5000),'r-','Label','Forelimb Landing');
    hold off;
    title({'Plantaris'; 'choose inactive first';'then active range'})
    ylabel('Activity (V)')
    xlabel('frame')
    legend('Raw','Smoothed','Location','northeast')
    % choose 4 points for range when plantaris & anconeus are at rest
    [Pl, ~] = ginput(4);
    close(figName);
    start_plant_inact = round(Pl(1));
    last_plant_inact = round(Pl(2));
    start_plant_act = round(Pl(3));
    if ((round(Pl(4))-5000)/10000) > times{1,3} % if plantaris activity extends beyond hindlimb liftoff
        last_plant_act = times{1,3}*10000 + 5000; % set cutoff at hindlimb liftoff
    else % if turns off before HL liftoff then keep as input
        last_plant_act = round(Pl(4));
    end

    figName = [R ': Select ranges of Anconeus, inactive then active'];
    figure('Name',figName)
    plot(E.anconeus,'k-')
    hold on;
    plot(EMG_smooth(:,2),'r-')
    xline(5001,'k--','Label','Start Moving')
    xline((times{1,2}*10000 + 5000),'Label','Forelimb Liftoff');
    xline((times{1,3}*10000 + 5000),'b-','Label','Hindlimb Liftoff');
    xline((times{1,4}*10000 + 5000),'r-','Label','Forelimb Landing');
    hold off;
    title({'Anconeus'; 'choose inactive first';'then active range'})
    ylabel('Activity (V)')
    xlabel('frame')
    % choose 4 points for range when plantaris & anconeus are at rest
    [An, ~] = ginput(4);
    close(figName);
    start_anc_inact = round(An(1));
    last_anc_inact = round(An(2));
    start_anc_act = round(An(3));
    if ((round(An(4))-5000)/10000) > times{1,4} % if plantaris activity extends beyond hindlimb liftoff
        last_anc_act = times{1,4}*10000 + 5000; % set cutoff at touchdown
    else % if turns off before touchdown then keep as input
        last_anc_act = round(An(4));
    end

    % Saving start and stop times when muscle is inactive in table
    start_inactive = [start_plant_inact; start_anc_inact];
    last_inactive = [last_plant_inact; last_anc_inact];
    start_active = [start_plant_act; start_anc_act];
    last_active = [last_plant_act; last_anc_act];
    t_muscles = table(start_inactive,last_inactive,start_active,last_active,'RowNames',{'Plantaris','Anconeus'});

else % have activation timing of muscles so use specified values
    start_plant_inact = round(varargin{1,1}.start_inactive('Plantaris'),0); % rounding ensures integer & supression of warning if number is in scientific format
    last_plant_inact = round(varargin{1,1}.last_inactive('Plantaris'),0);
    start_plant_act = round(varargin{1,1}.start_active('Plantaris'),0);
    last_plant_act = round(varargin{1,1}.last_active('Plantaris'),0);
    start_anc_inact = round(varargin{1,1}.start_inactive('Anconeus'),0);
    last_anc_inact = round(varargin{1,1}.last_inactive('Anconeus'),0);
    start_anc_act = round(varargin{1,1}.start_active('Anconeus'),0);
    last_anc_act = round(varargin{1,1}.last_active('Anconeus'),0);
    t_muscles = varargin{1,1};
end

%% Plantaris
% Determining Plantaris average resting activity and SD of raw data (should be
% >=2 SD's of rectified)
P_rest_avg = mean(EMG_smooth(start_plant_inact:last_plant_inact,1));
P_sd = std(E{start_plant_inact:last_plant_inact,1});

% Defining Plantaris Threshold
Plant_thresh = P_rest_avg + 2*P_sd;

% Activation timing of Plantaris
Plant_Active = EMG_smooth(start_plant_act:last_plant_act,1) >= Plant_thresh;
Plant_ON_takeoff = (find(Plant_Active, 1, 'first') + start_plant_act - 1)/10000 - 0.5 - times{1,3}; % units of seconds
Plant_ON_start = (find(Plant_Active,1,'first') + start_plant_act - 1)/10000 - 0.5; % units of seconds
Plant_act_dur_to = (find(Plant_Active,1,'last') - find(Plant_Active,1,'first'))/10000; % duration of activity in seconds


%% Anconeus
% Determining Anconeus average resting activity and SD of raw data (should be
% >=2 SD's of rectified)
A_rest_avg = mean(EMG_smooth(start_anc_inact:last_anc_inact,2));
A_sd = std(E{start_anc_inact:last_anc_inact,2});

% Defining Anconeus Threshold
Anc_thresh = A_rest_avg + 2*A_sd;

% Activation of Anconeus
Anc_Active = EMG_smooth(start_anc_act:last_anc_act,2) >= Anc_thresh;
if start_anc_act >= times{1,3}*10000+5000 % when initial activation of anconeus was marked after takeoff
    idx_Anc_activates = find(EMG_smooth(start_anc_act:last_anc_act,2) >= Anc_thresh,1,"first") + (start_anc_act - 1); % when anconeus first activates following hindlimb liftoff (assuming inactive before HL liftoff)
    Anc_Active_air = EMG_smooth(start_anc_act:last_anc_act,2) >= Anc_thresh; % use that marker
else % initial activation of anconeus was marked prior to hindlimb takeoff & we only care about aerial phase & just prior to takeoff
    idx_last_inactive = find(EMG_smooth(1:round(times{1,3}*10000+5000,0),2) < Anc_thresh,1,"last"); % index where muscle is last inactive prior to takeoff
    idx_Anc_activates = idx_last_inactive + 1;
    Anc_Active_air = EMG_smooth(round(times{1,3}*10000+5000,0):last_anc_act,2) >= Anc_thresh;
end
Anc_ON_HLliftOff = (idx_Anc_activates)/10000 - 0.5 - times{1,3}; % time in seconds Anconeus turns ON relative to Hindlimb Takeoff
Anc_act_dur_air = (find(Anc_Active_air,1,'last') - find(Anc_Active_air,1,'first'))/10000; % duration of activity in seconds

end

%% Trim_Igor (function)
function E = Trim_Igor(Igor,times)
% This function trims Igor files down to time of interest, 60s after
% trigger's falling edge

% INPUTS
% Igor - table of data collected from Igor containing plantaris activity,
%        anconeus activity and the trigger in that order for columns.
% times - table of times where first row is relative to within jump and
%         second row is relative to entire 60s recording.

% OUTPUTS
% E - table containing the raw EMG values for each muscle

I = table2array(Igor);
Triggered = find(I(:,3) < 1); % find index when trigger is depressed
% finding 60 second time frame of interest
I = I(Triggered(1,1):(Triggered(1,1)+599999),:); % only keeps 60s following falling edge of trigger

% Start time will be 0.5 seconds before onset of movement
start = round(times{2,1}*10000 + 1) - 5000; % times are in movie frames (250 fps) and need to be at 10,000Hz for EMG, subtract off 5000 frames to get 0.5s prior to movement initiation
stop = round(times{2,5}*10000 + 1);
if start <= 0 % if jump is early enough where there isn't 0.5s before start of recording
    Z_start = zeros((-1*start+1),1); % add zeros for frames that pre-exist recordign
    plantaris = [Z_start; I(1,stop,1)];
    anconeus = [Z_start; I(1,stop,2)];
else
    plantaris = I(start:stop,1); % raw voltage of plantaris
    anconeus = I(start:stop,2); % raw voltage of anconeus
end

E = table(plantaris,anconeus);
end