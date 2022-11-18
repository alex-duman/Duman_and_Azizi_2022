# Duman and Azizi (Submitted 2022)
Code repository for Duman and Azizi (submitted 2022) manuscript in the Journal of Experimental Biology.

## Table of Contents
* [Installation](##Installation)
* [Workflow](##-Workflow)
* [Code Repository (GitHub)](##Code-Repository-(GitHub))
* [Data Deposition (Dryad)](##Data-Deposition-(Dryad))

## Installation

Clone the repository to your device in order to use and reproduce the results from Duman and Azizi (submitted 2022). You should then be able to view files with Microsoft Excel and MATLAB. 

If you want to regenerate the csv files containing the processed data you will need to download the raw data files from Dryad (see Data Deposition section below). You do not need to download the folder containing all video recordings to regenerate these files, but they are provided for your reference. Make sure the downloaded data folders (Digitized Files, Figures, Igor Files, and Reflex Data) are within the same directory that you cloned the directory to.

[(Back to Top)](#Duman-and-Azizi-(Submitted-2022))

## Workflow
1. Start by running **Analyze_ALL_Reflex_Trials.m** to reproduce *Analyzed_Data_Reflex.csv* and get results of muscle stretch reflex experiment.
2. Then run **Statistics_Reflex_FINAL.m** to reproduce plots from *Figure 3* & Table_Reflex which can be compared to Sheet 1 of Compiled_Statistics.xlsx
3. Next, run **Analyze_Script_FINAL.m** which should allow you to save results from the continuous jumping experiment as *Analyzed_Data.csv* and should match the provided *Analyzed_Data_Final.csv*.
4. Then run **Statistics_FINAL.m** using the provided Analyzed_Data_Final.csv or your Analyzed_Data.csv (be sure to change the appropriate line where D_all data is read in to run your version). This replicates *Figures 4-6* and *Table 1* in the manuscript. Variables Table_1 and Table_EMG can also be compared to Sheet 1 and Sheet 2 of Compiled_Statistics.xlsx.
5. Next, run **Seqential_Bonferroni.m** which uses Compiled_Statistics.xlsx to then create a ranking to determine the significance of the multiple statistical tests performed throughout the manuscript. To view the results check the variable *Seq_Bon_ranked*.
6. Finally, to reproduce *Supplemental Figure 1* run the **Correlating_JumpDist_w_ElbExt_td.m** script.

[(Back to Top)](#Duman-and-Azizi-(Submitted-2022))

## Code Repository (GitHub)
Below is a brief description of each file contained within the repository. More information including descriptions of specific inputs and outputs are provided within the files themselves.

**Analysis_and_Plot_Reflexes.m** - function that reads in raw electromyography (EMG) and kinematic data recorded from muscle stretch reflex trials, processes it and produces figures if Plot term is set to 'Y' or 'Yes' (otherwise graphs will not be output), it will always output a table containing variables of interst.

**Analyze_ALL_Reflex_Trials.m** - this file analyzes all muscle stretch reflex trials by iterating through the list of trials in List_of_Reflex_Trials.xlsx and passing each trial into the Analysis_and_Plot_Reflexes. It then compiles data from all muscle stretch reflex trials into a table of results, and saves this table to Analyzed_Data_Reflex.csv.

**Analyze_Script_FINAL.m** - This file analyzes all the jump trials recorded in the Hindlimb_Proprioception_Timing.xlsx spreadsheet. It will also load the t_muscle.mat file which has information about the activation timing of the anconeus (elbow extensor) and plantaris (ankle extensor) muscles. If you would prefer to select your own onset and offest times for these muscles delete the t_muscle.mat file and the code will walk you through selecting the appropriate times by visualizing the data for each trial. This code will also prompt you to save and overwrite the resulting table Analyzed_Data_Final.csv and t_muscle.mat file.

**Analyzed_Data_Final.csv** - a table with the variables of interest from the continuous jumping experiments derived from analyses in Analyze_Script_FINAL.m (should have same structure as Analyze_Data.csv output).
* **Date** = date and corresponding video folder associated with when the jump was recorded and where the corresponding video is stored (column 1; see Data Deposition section for more about video data)
* **trial** = the trial name corresponding to the video file name the jump was collected from (column 2)
* **condition** = the experimental condition where pre is pre-surgical intervention, sham is 1 week following the sham surgery where the sciatic nerves were left intact, and post is 6 months post-reinnervation surgery (column 3)
* **Jump** = the unique number identifying each jump within a given trial, ordered sequentially from when the occured within the video trial (column 4)
* **toad** = the unique number assigned to individual toads (*Rhinella marina*; column 5)
* **avg_dElbExt** = average rate of elbow extension in units of degrees per second (column 6)
* **ElbExt_td** = elbow extension at time of forelimb touchdown in units of degrees (column 7)
* **time_ElbExt_to** = time elbow extension begins following hindlimb takeoff in units of seconds (column 8)
* **time_ElbExt_td** = time elbow extension begins prior to forelimb touchdown in units of seconds (column 9)
* **avg_dhLER_to** = average rate of hindlimb extension (ratio) in units of limb lengths per second (column 10; here a positive rate implies limb extension and negative rate is limb flexion)
* **max_dhLER_to** = maximal rate of hindlimb extension during takeoff in units of limb lengths per second (column 11)
* **dhLER_to** = the instantaneous rate of hindlimb extension at the moment of hindlimb takeoff in units of limb lengths per second (column 12)
* **t_HLflexes** = time the hindlimb begins flexing during the aerial phase relative to hindlimb takeoff in units of seconds (column 13)
* **avg_belt_vel** = average treadmill belt velocity (or speed, positive velocity is treadmill moving backward so animal is jumping forward; column 14)
* **takeoff_velocity** = the instantaneous speed of the center of mass at time of hindlimb takeoff in units of meters per seconds (column 15)
* **jump_dist** = the distance covered during the jump in units of meters (column 16)
* **Plant_act_dur_to** = the duration of *plantaris* muscle (ankle extensor) activation prior to hindlimb takeoff in units of seconds (column 17)
* **Anc_ON_HLliftOff** = the duration of time between hindlimb liftoff and when the *anconeus* muscle (elbow extensor) became active in units of seconds (negative values represent trials where the muscle became active prior to hindlimb liftoff and positive values when the muscle became active after hindlimb liftoff; column 18)
* **Anc_act_dur_air** = the duration of *anconeus* muscle activation during the aerial phase of the jump in units of seconds (column 19)
* **Crash** = represents whether the landing behavior was characterized by coordinated forelimb-mediated landing (crash = 0) or if the toad allowed its torso and or head to crash into the ground (crash = 1; column 20)

**Analyzed_Data_Reflex.csv** - a table with the variables of interest from the muscle stretch reflex experiments derived from analyses in Analyze_ALL_Reflex_Trials.m.
* **Individual** = corresponds to the number associated with each individual toad (column 1)
* **Condition** = the experimental condition where Pre-op is pre-surgical intervention, Sham is 1 week following the sham surgery where the sciatic nerves were left intact, 1wk Post-op is 1 week following nerve reinnervation surgery, 3mo Post-op is 3 months following nerve reinnervation, 6mo Post-op is 6 months following nerve reinnervation, and Post-mortem is recordings taken after the animals have passed (column 2)
* **Trial** = corresponds to the number associated with replicating the muscle stretch reflex measurement, with three replicate recordings being taken for each timepoint and eventually averaged in later analyses (remove pseudoreplication; column 3)
* **EMG_thresh** = the threshold for determining activation onset from the *plantaris* muscle (ankle extensor) activity in units of volts (see Analysis_and_Plot_Reflexes.m for more detail on how EMG_thresh is calculated; column 4)
* **I_flexion** = the electromyographic (EMG) intensity (area under the rectified EMG signal) during the period of ankle flexion in units of Volt seconds (column 5)
* **I_extension** = the EMG intensity (area under the rectified EMG signal) during the period of ankle extension in units of Volt seconds (column 6)
* **Ratio_I_ext_flex** = the ratio of I_extension over I_flexion which is a unitless way of determining the relative amount of activation during ankle extension relative to ankle flexion; in a muscle with intact stretch reflex we would expect this ratio to be much greater than 1 (one; column 7)
* **t_on_flex** = the duration of time during the period of ankle flexion that the *plantaris* muscle is active in units of seconds (column 8)
* **t_on_ext** = the duration of time during ankle extension the muscle is active in units of seconds (column 9)
* **Ratio_t_on_ext_flex** = the ratio of t_on_ext over t_on_flex which is a unitless way of determining the relative amount of activation during ankle extension relative to ankle flexion; in a muscle with intact stretch reflex we would expect this ratio to be much greater than 1 (one; column 10)

**Compiled_Statistics.xlsx** - contains the statistical tables from the Statistics_FINAL.m with Sheet 1 relating to Table_1 regarding forelimb and hindlimb behavior (Table 1 in the manuscript) and Sheet 2 representing the Table_EMG relating to the electromyography results.
* **Sheet 1** = results from Linear Mixed Effects (LME) Models obtained from analyzing hypotheses regarding muscle stretch reflex as well as the forelimb and hindlimb behavior
    * **Variable** = the dependent variable of interest (column 1)
    * **Condition** = the corresponding experimental conditions within each model (column 2)
    * **Mean** = the group mean (average) for each condition (column 3)
    * **95% CI** = the 95% Confidence Interval of the mean, we are 95% confident that the group mean exists within the reported 95% CI value away from our observed mean (column 4)
    * **SD** = the standard deviation for each group (column 5)
    * **DF1** = the numerator degrees of freedom for the dependent variable of interest (column 6)
    * **DF2** = the denominator degrees of freedom for the dependent variable of interest (column 7)
    * **F_stat** = the corresponding F-statistic representing the variation between groups relative to the variation within groups (column 8)
    * **p-value** = the corresponding p-value derived from the F-statistic (column 9)
    * *NOTE: some rows contain empty parentheses because these values would be redundant, for example rows 2-7 all correspond to the (Ratio)_I_ext_flex variable*
* **Sheet 2** = results from the Student's T-Tests for comparing experimental results of electromyography (muscle acitivty) 6 months following nerve reinnervation with values obtained from intact toads reported in pre-existing literature
    * **Variable** = the dependent variable of interest (column 1)
    * **Mean_post** = the mean (average) value recorded from the experimentally derived 6 month post-reinnervation toads (column 2)
    * **95%_CI_post** = the 95% Confidence Interval of the mean for individuals 6 months following reinnervation surgery, we are 95% confident that the group mean exists within the reported 95% CI value away from our observed mean (column 3)
    * **Mean_Ref1** = the mean (average) value reported from the first reference paper (column 4)
    * **95%_CI_Ref1** = the 95% Confidence Interval of the mean for the first reference paper (if sufficient data was provided; column 5)
    * **p_val_Ref1** = the p-value obtained from a One-Sample Student's T-Test comparing the 6 month post-reinnervated data with the reported mean from the first reference paper (column 6)
    * **Mean_Ref2** = the mean (average) value reported from the second reference paper if applicable (column 7)
    * **95%_CI_Ref2** = the 95% Confidence Interval of the mean for the second reference paper (if applicable and sufficient data was provided; column 8)
    * **p_val_Ref2** = (may be typo of p_val_Ref3) the p-value obtained from a One-Sample Student's T-Test comparing the 6 month post-reinnervated data with the reported mean from the second reference paper (column 9)

**Correlating_JumpDist_w_ElbExt_td.m** - script to determine whether elbow angle (extension) at touchdown correlates with jump distance. The visualization produced reflects Supplementary Figure 1 associated with the manuscript.

**Hindlimb_Proprioception_Timing.xlsx** - a table with information pertaining to the relative timing of specific components of each jump and whether or not the jump resulting in a coordinated (0) or crash (1) landing. Timing data has units of video frames used to identify key times (video data collected prior to September 2020 was recorded at 100 fps and videos from September 2020 were recorded at 250 fps).
* **Date** = date and corresponding video folder associated with when the jump was recorded and where the corresponding video is stored (column 1; see Data Deposition section for more about video data)
* **trial** = the trial name corresponding to the video file name the jump was collected from (column 2)
* **Jump** = the unique number identifying each jump within a given trial, ordered sequentially from when the occured within the video trial (column 3)
* **Start** = the video frame number corresponding to when toad started moving (movement initiation; column 4)
* **Forelimb Liftoff** = frame number corresponding to when the toad's forelimb was last in contact with the ground (column 5)
* **Hindlimb Liftoff** = frame number corresponding to when the toad's hindlimb were last in contact with the ground (column 6)
* **Forelimb Touchdown** = frame number corresponding to when the toad's forelimb first contacted the ground to initiate landing (column 7)
* **Body Touchdown** = frame number corresponding to when any body part other than the forelimbs made contact with the ground and marks the end of forelimb-only-mediated landing (column 8)
* **Crash (Y=1, N=0)** = represents whether the landing behavior was characterized by coordinated forelimb-mediated landing (0, or No crash) or if the toad allowed its torso and or head to crash into the ground (1 or Yes crashed; column 9)
* *NOTE: trials recorded prior to September 2020 were recorded at 100 fps and trials recorded during September 2020 were recorded at 250 fps.*

**List_of_Reflex_Trials.xlsx** - lists trial information for the muscle reflex experiments.
* **File** = the current file name for the corresponding muscle stretch reflex trial (column 1)
* **Iniital Trial Name** = original file names associated with video and Igor recordings with values differing from the current File naming convention highlighted in red (column 2)
* **date** = the date when the muscle stretch reflex trial was performed

**Sequential_Bonferroni.m** - this file uses the Compiled_Statistics.xlsx workbook to perform a Sequential Bonferroni ranking in order to account for the multiple tests throughout the manuscript. The Seq_Bon_ranked is the variable that lists each dependent variable (1 for each hypothesis test), corresponding p-value, the threshold for rejection, and whether or not we can claim a significant difference given our multiple hypothesis tests.

**Statistics_FINAL.m** - this uses the provided Analyzed_Data_Final.csv to calculate the results of the Linear Mixed Effects (LME) Models and One-Sample Student's T-tests to replicate Figures 4-6 in the manuscript and Table 1. To perform the Sequential Bonferroni Table_1 and Table_EMG are coppied into Sheet 1 and Sheet 2 of Compiled_Statistics.xlsx, respectively (results from Statistics_Reflex_FINAL.m are also included in Sheet 1). 

**Statistics_Reflex_FINAL.m** - uses Analyzed_Data_Reflex.csv to reproduce plots from Figure 3 & Table_Reflex which is included at the top of Sheet 1 in Compiled_Statistics.xlsx for consideration with the Sequential Bonferroni.

**t_muscle.mat** - mat file that contains the variable t_muscle which is a column-cell array where each row represents a different jump trial. Within each cell corresponding to a trial where muscle activity was recorded (trials in September 2020, 6 months post-reinnervation surgery) there exists a table with activation onset and offset for the two muscles of interest. This file is optional and you can delete it to try assigning your own onset and offset times for each muscle, however these values are necessary to replicate the results presented in the manuscript.




## Data Deposition (Dryad)
All the necessary data files to perform analyses and produce the results, figures and tables from the manuscipt can be found in this deposition. To access the code necessary to generate the results visit: https://github.com/alex-duman/Duman_and_Azizi_2022

[(Back to Top)](#Duman-and-Azizi-(Submitted-2022))

### Folders, Data Files & Descriptions
**Digitized Files** - contains the 3D coordinates of the toad's limbs for each continuous hopping trial. There are two files corresponding to each trial, the first having the trial name followed by '_xyzpts' and the second having '_TH_xyzpts' as the suffix to the trial name (e.g. Rmar01_01_xyzpts.csv & Rmar01_01_TH_xyzpts.csv). Data in this folder are used to perform statistical analyses and create the figures, tables and results for the manuscript.

* **TrialName_xyzpts.csv** - (ex: Rmar01_xyzpts.csv; only applies to trials prior to 6 months post-reinnervation; prior to September 2020) contain the x, y and z coordinates of the body across columns with each row representing time or a unique video frame through the continuous hopping trial containing only jumps of interest. These files relate to the videos with corresponding trial names within the Videos > Continuous Hopping >... subfolders corresponding to the time points prior to 6 months post-reinnervation. This data was derived from videos recorded with a frame rate of 250 fps using the Hedrick Lab's DLTdv digitizing tool (https://biomech.web.unc.edu/dltdv/). 
    * **pt1** = a point on the treadmill belt that was tracked in units of meters (_X, _Y, _Z in columns 1-3)
    * **pt2** = position of animal's right wrist in meters (_X, _Y, _Z in columns 4-6)
    * **pt3** = position of animal's right elbow in meters (_X, _Y, _Z in columns 7-9)
    * **pt4** = position of animal's right shoulder in meters (_X, _Y, _Z in columns 10-12)
    * **pt5** = position of the animal's longest toe tip on right hindlimb in meters (_X, _Y, _Z in columns 13-15)
    * **pt6** = position of animal's right metatarsophalangeal joint in meters (_X, _Y, _Z in columns 16-18)
    * **pt7** = position of animal's right ankle in meters (_X, _Y, _Z in columns 19-21)
    * **pt8** = position of animal's right knee in meters (_X, _Y, _Z in columns 22-24)
    * **pt9** = position of animal's posterior protrusion of the ischium in meters (_X, _Y, _Z in columns 25-26)
* **TrialName_xyzpts.csv** - (ex: Rmar01_01_xyzpts.csv; only applies to trials during 6 months post-reinnervation; September 2020) contain the x, y and z coordinates of the body across columns with each row representing time or a unique video frame through the continuous hopping trial containing only jumps of interest. These files relate to the videos with corresponding trial names within the Videos > Continuous Hopping >... subfolders corresponding to time points at 6 months post-reinnervation. This data was derived from videos recorded with a frame rate of 250 fps using the Hedrick Lab's DLTdv digitizing tool (https://biomech.web.unc.edu/dltdv/). 
    * **pt1** = a point on the treadmill belt that was tracked in units of meters (_X, _Y, _Z in columns 1-3)
    * **pt2** = position of animal's right wrist in meters (_X, _Y, _Z in columns 4-6)
    * **pt3** = position of animal's right elbow in meters (_X, _Y, _Z in columns 7-9)
    * **pt4** = position of animal's right shoulder in meters (_X, _Y, _Z in columns 10-12)
    * **pt5** = position of animal's right metatarsophalangeal joint in meters (_X, _Y, _Z in columns 13-15)
    * **pt6** = position of animal's right ankle in meters (_X, _Y, _Z in columns 16-18)
    * **pt7** = position of animal's right knee in meters (_X, _Y, _Z in columns 19-21)
        * **toe & ischium** = NOT included in this file but in a separate file with name convention like 'TrialName_TH_xyzpts.csv' (see below)
* **TrialName_TH_xyzpts.csv** - (ex: Rmar01_01_TH_xyzpts.csv) contain the x, y and z coordinates of the toe and hip (estimated as the posterior protrusion of the ischium) in columns with rows representing time or a unique video frame through the continuous hopping trial containing only jumps of interest. This file only exists for trials recorded during September 2020 and therefore only exists for 6 month post-reinnervation trials (toe and hip data are included in the main file for trials recorded prior to 6 months post-reinnervation). This data was derived from videos recorded with a frame rate of 250 fps (all trials during 6 months post-reinnervation time point) using the Hedrick Lab's DLTdv digitizing tool (https://biomech.web.unc.edu/dltdv/).
    * **pt1** = position of the animal's longest toe tip on right hindlimb in meters (_X, _Y, _Z in columns 1-3)
    * **pt2** = position of animal's posterior protrusion of the ischium in meters (_X, _Y, _Z in columns 4-6)

**Igor Files** - contains the electromyography (EMG) data recorded in Igor for trials recorded 6 months post-reinnervation (e.g. Rmar01_01_emg.csv).
* **TrialName_emg.csv** - (ex: Rmar01_01_emg.csv) contain the electromyographic (EMG) and other data through time collected during the 6 months post-reinnervation continuous hopping (trials associated with videos in the Videos > Continuous Hopping > Post-Reinnervation... subfolders). This data was recorded for 65 seconds and had a sampling frequency of 10 kHz.
    * **plantaris** = muscle activity (EMG) in Volts recorded from plantaris muscle (ankle extensor, column 1)
    * **anconeus** = muscle activity (EMG) in Volts recorded from the anconeus muscle (elbow extensor, column 2)
    * **trigger** = recording of the cameras' trigger signal in units of Volts which remains at approximately +3V until trigger drops signal to approximately 0V, and triggers the start of recording on the falling edge (column 3, video recording only occurs for 60 seconds at either 100 fps or 250 fps)

**Reflex Data** - contains the Igor data (e.g. Rmar01_01_Igor.csv) and 2D coordinates of the ankle joint that allow for the analysis of the muscle stretch reflex data.
* **TrialName_Igor.csv** - (ex: Rmar01_01_Igor.csv) contain the electromyography (EMG), motor force, motor length, and trigger through time to synchronize the Igor data and video corresponding the muscle stretch reflex test. This data was recorded for 5 seconds and had a sampling frequency of 10 kHz.
    * **EMG** = muscle activity (electromyography) in Volts of the plantaris muscle (ankle extensor) through time (column 1)
    * **force** = the resistive force in Newtons applied and measured by the motor (used to accelerate the foot to enduce rapid ankle flexion, column 2)
    * **length** = relative position of the motor in units of Volts (column 3)
    * **trigger** = recording of the cameras' trigger signal in units of Volts which remains at approximately +3V until trigger drops signal to approximately 0V, and triggers the camera to save the previous 1 second of video recorded (column 4)
* **TrialName_xypts.csv** - (ex: Rmar_01_01_xypts.csv) contain the 2D (x and y) coordinates from the muscle stretch reflex test videos. Columns represent the x and y positions of points tracked (listed below) and rows represent time or a unique video frame through the duration of the reflex trial. This data was derived from videos recorded for 1 second with a frame rate of 1964 fps using the Hedrick Lab's DLTdv digitizing tool (https://biomech.web.unc.edu/dltdv/).
    * **pt1_cam1** = left-most corner of table that acted as stationary reference (_X and _Y, columns 1-2)
    * **pt2_cam1** = right-most corner of table that acted as stationary reference, together with pt1_cam1 could form a vector to quantify the angular orientation of the shank (strapped to table) and ultimately track the change in the ankle joint (_X and _Y columns 3-4)
    * **pt3_cam1** = proximal corner on hinge attached to foot (_X and _Y, columns 5-6)
    * **pt4_cam1** = distal corner on hinge attached to foot, together with pt3_cam1 could form a vector to quantify the angular orientation of the foot (strapped to hinge) and ultimately track the change in the ankle joint (_X and _Y columns 7-8)

**Videos** - contains all the videos from this project separated into the portion devoted to measuring the Muscle Stretch Reflex & Continuous Hopping experiments.

**Videos > Continuous Hopping** - contains subfolders that hold videos (.avi files) for each time toads were recorded hopping (e.g. Pre-Surgery (20_Jan_2020A), Post-Reinnervation (25_Sep_2020), etc.). All videos in folders dated prior to 25th September 2020 were filmed at 100 fps, and videos recorded on or after 25th of September 2020 were recorded at 250 fps. Associated calibration files (csv files) are provided and correspond to the Hedrick Lab's DLTdv digitizing tool (https://biomech.web.unc.edu/dltdv/).
* **Cal_Date_DLTcoefs.csv** - (ex: Post-Reinnervation(25_Sep_2020) > Cal_25_Sept_20_DLTcoefs.csv) contain the coefficients necessary to calibrate the two camera views (front and side) for a trial in order to perform the direct linear transform (DLT) and compute 3D coordinates in cartesian (x, y, z) space. It is important to note that we calibrateed all videos by loading the front video before the side video. For more detail visit (https://biomech.web.unc.edu/dltdv/).
* **TrialName_front.avi** - (ex: Post-Reinnervation(25_Sep_2020) > Rmar01_01_front.avi) the video taken from the camera positioned toward the front (slightly frontal plane view) of the treadmill during the continuous hopping trial. Videos were recorded at 100 fps prior to September 2020 and at 250 fps during September 2020 (6 months post-reinnervation trials). Original videos were 60 seconds in length however these videos have been trimmed to only include jumps were the animal did not jump directly into side walls in the order of occurance.
* **TrialName_side.avi** - (ex: Post-Reinnervation(25_Sep_2020) > Rmar01_01_side.avi) the video taken from the camera positioned on the side (lateral view) of the treadmill during the continuous hopping trial. Videos were recorded at 100 fps prior to September 2020 and at 250 fps during September 2020 (6 months post-reinnervation trials). Original videos were 60 seconds in length however these videos have been trimmed to only include jumps where the animal did not jump directly into side walls and remove any periods where the animal was not jumping. Jumps appear in the order of occurance during each trial's video.

**Videos > Muscle Stretch Reflex** - contains subfolders that hold videos (.avi files) for each time the muscle stretch reflex was measured (e.g. Pre-surgery, Sham-surgery (1 week), Post-mortem, etc.). All of these videos were filmed at 1964 fps.
* **TrialName.avi** - (ex: Post-mortem > Rmar_01_d1.avi) the video recording of the muscle stretch reflex test on the plantaris muscle (ankle extensor). Videos were recorded for 1 second at a frame rate of 1964 fps.

[(Back to Top)](#Duman-and-Azizi-(Submitted-2022))
