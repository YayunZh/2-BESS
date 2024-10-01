# LS-HiPPP-EnergyDesign

I. To run the code, you just need to open the main.m file and run it section by section. The utimate goal is to maximize the energy utilization of battery networks. The main.m is divided into serveral sections: User Input, Task 4 Energy Design Layer 1 Converters, Task 4-1: Energy Design Layer 2 Converters, Task 4-4: Compare to C-PPP, Task 4-5: Compare to C-PPP and FPP
***************************************************************************************************************************************************
II. Tunable parameters in each sections
A. User Input 
votlage of ith batteries: Input{i}.volt
charge of ith batteries: Input{i}.qlim
mean of the charge of ith batteries: Input{i}.qlim_mu
variance of the charge of ith batteries: Input{i}.qlim_var
number of layer 1 convertetrs: ENV.Avg_Conv.Num
number of layer 2 convertetrs: ENV.Var_Conv.Num
adjunct matrix of the direct power transfer link: P_direct_mat_in
adjunct matrix of the differential power link: P_diff_var_mat_in

B. Task 4 Energy Design Layer 1 Converters
the numer of random searches to find the optimal placement of averaging converter: ENV4.Avg_Conv.trial_num (typical number is 10000, to save time can decrease to 1000)
how to partition the layer 1 converters to achieve the economic of scales: ENV4.Avg_Conv.partition (typical number is 1)

C. Task 4-1: Energy Design Layer 2 Converters 
number of battery distributions to sweep: ENV4_1.Sweep.Stat.Bat, by assuming battery follows Gaussian, the following two vectors governs teh battery distributions
mean vector of battery to sweep: ENV4_1.Sweep.Bat{i}.qlim_mu
std vector of battery to sweep: ENV4_1.Sweep.Bat{i}.qlim_var
number of converter ratings to sweep: ENV4_1.Sweep.Stat.Conv
number of trials for MC simulation: (typical number is 100, to save time can decrease to 20)

D. Task 4-4: Compare to C-PPP
number of trials for MC simulation: ENV4_trad.Var_Conv.MC_trial (typical number is 100, to save time can decrease to 20)

E. Task 4-5: Compare to C-PPP and FPP
number of trials for MC simulation: ENV4_fpp.Var_Conv.MC_trial (typical number is 100, to save time can decrease to 20)

****************************************************************************************************************************************************
III. Guide to run each task
A. Guide for Task 4-1: Run User input -> Run Task 4 Energy Design Layer 1 Converters -> Run Task 4-1: Energy Design Layer 2 Converters -> Run Plot Task 4-1; 

B. Guide for Task 4-4: Run User input -> Run Task 4 Energy Design Layer 1 Converters -> Run Task 4-1: Energy Design Layer 2 Converters -> Run Task 4-4: Compare to C-PPP -> Run Plot Task 4-4; 

C. Guide for Task 4-4: Run User input -> Run Task 4 Energy Design Layer 1 Converters -> Run Task 4-1: Energy Design Layer 2 Converters -> Run Task 4-4: Compare to C-PPP -> Task 4-5: Compare to C-PPP and FPP -> Run Plot Task 4-5
