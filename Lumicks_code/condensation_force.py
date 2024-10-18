"""
Script to get the topo condensation force.
"""
import numpy as np
import lumicks.pylake as lk
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
%matplotlib qt 

filepath = r"Y:\Lumicks\C-trap user data\Meiling\240628_121032_unlabeled_HsTOP2A_2p5nM_No_Sytox_7_5per_No_flow_load_hold 6um"
filename = "unlabeled_HsTOP2A_2p5nM_No_Sytox_7_5per_No_flow_load_hold 6um tether 14 Topo Kymo.h5"
        
file = lk.File(filepath + "\\" + filename)
        
fig_title = ("Condensation force analysis")

# The condensation force was measured during a kymo in the Lumicks data.  We will use the Kymo to extract the data.
for key, value in file.kymos.items():
    keyVal = key;
kymo = file.kymos[keyVal];

# The following conversion is required to get high frequency data from within the kymograph
# This alignment removes all data from before the kymograph starts from the HF/LF data
kymo_timestamp_ranges = kymo.line_timestamp_ranges()
kymo_start_timestamp = kymo_timestamp_ranges[0][0] # This is the first timepoint with kymo data (approximately)
kymo_end_timestamp = kymo_timestamp_ranges[-1][-1]

# Data from the entire marker
# Marker has 1 kymographs
distance_raw = file["Distance"]["Distance 1"].data
distance_time_raw = file["Distance"]["Distance 1"].seconds

force_raw = file["Force LF"]["Force 2x"].data
forceTime_raw = file["Force LF"]["Force 2x"].seconds # Time in seconds

force_distance_down = file["Force HF"]["Force 2x"].downsampled_to(20, method="force").data
dist_distance_down = file["Distance"]["Distance 1"].downsampled_to(20, method="force").data
distance_time_downsampled = file["Distance"]["Distance 1"].downsampled_to(20, method="force").seconds
marker_timestamps = file["Distance"]["Distance 1"].downsampled_to(20, method="force").timestamps # Dist and Force are taken at same timestamps
marker_timestamps_raw = file["Distance"]["Distance 1"].timestamps

if len(dist_distance_down) > len(force_distance_down):
    while len(dist_distance_down) > len(force_distance_down):
        dist_distance_down = dist_distance_down[:-1]
        distance_time_downsampled = distance_time_downsampled[:-1]
        
if len(dist_distance_down) < len(force_distance_down):
    while len(dist_distance_down) < len(force_distance_down):
        force_distance_down = force_distance_down[:-1]

kymo_start_loc_down = np.where(marker_timestamps > kymo_start_timestamp)[0][0]
kymo_start_loc_raw = np.where(marker_timestamps_raw > kymo_start_timestamp)[0][0]
kymo_start_time_down = distance_time_downsampled[kymo_start_loc_down]

kymo_end_loc_down = np.where(marker_timestamps < kymo_end_timestamp)[0][-1]
kymo_end_time_down = distance_time_downsampled[kymo_end_loc_down]

hold_start_time = kymo_start_time_down
hold_end_time = hold_start_time + 99


hold_slice=np.where((distance_time_downsampled > hold_start_time)&(distance_time_downsampled<hold_end_time))
   
###############################################################################
#Force generated during hold, get a mean of the force of last 10s (hold_force)
###############################################################################
hold_slice_last_10s = np.where((distance_time_downsampled>(hold_end_time-10)) &(distance_time_downsampled<hold_end_time))
hold_force = np.mean(force_distance_down[hold_slice_last_10s])


time_cor_down = distance_time_downsampled[kymo_start_loc_down]
time_cor_raw = distance_time_raw[kymo_start_loc_raw]


# This code starts the marker immediately before starting the kymo
# We can correct for this slightly
distance_time_downsampled = distance_time_downsampled - time_cor_down
dist_down_hold = dist_distance_down[distance_time_downsampled < hold_end_time]
time_down_hold = distance_time_downsampled[distance_time_downsampled < hold_end_time]
kymo_start_down_cor = kymo_start_time_down - time_cor_down
kymo_end_down_cor = kymo_end_time_down - time_cor_down    

distance_time_raw = distance_time_raw - time_cor_raw
forceTime_raw = forceTime_raw - time_cor_raw
distance_raw_hold = distance_raw[distance_time_raw < hold_end_time]
time_raw_hold = distance_time_raw[distance_time_raw < hold_end_time]
force_raw_hold = force_raw[forceTime_raw < hold_end_time]
hold_start_time = hold_start_time- time_cor_raw
hold_end_time = hold_end_time - time_cor_raw

# Make plot of force and extension to visualize trace
fig = plt.figure(layout="constrained")
gs1 = grd.GridSpec(2, 1, figure=fig)

axExt = fig.add_subplot(gs1[0])
axForce = fig.add_subplot(gs1[1])

fig.set_figheight(14)
fig.set_figwidth(18)

fig.suptitle(filename)
plt.subplots_adjust(top=0.925)
fig_title = fig_title

const_ylim_ext = [13, 4]
const_ylim_force = [-.5, 2]
time_xlim = [kymo_start_down_cor, kymo_end_down_cor]                     
    
# Make the extension plot                 
axExt.plot(distance_time_downsampled[hold_slice],dist_distance_down[hold_slice],'-', color = 'blue')
axExt.plot(distance_time_raw, distance_raw, '-', color = 'grey', alpha = 0.5)
axExt.set_title(f'Extension Plot (20 Hz) (Topo Channel)')
axExt.set_ylim(const_ylim_ext)
axExt.set_xlim(time_xlim)
axExt.set_ylabel("Extension (um)", fontsize =16)
axExt.set_xlabel("Time (s)",fontsize =16)

# Make the force plot
axForce.plot(distance_time_downsampled[hold_slice],force_distance_down[hold_slice],'-', color = 'blue')
axForce.plot(distance_time_downsampled[hold_slice_last_10s],force_distance_down[hold_slice_last_10s],'o',color = 'darkblue')
axForce.plot(forceTime_raw, force_raw, '-', color = 'grey', alpha = 0.5)
axForce.set_ylim(const_ylim_force)
axForce.set_xlim(time_xlim)
axForce.set_ylabel("Force (pN)",fontsize =16)
axForce.set_xlabel("Time (s)",fontsize =16)

axForce.text(np.max(distance_time_downsampled) + 0,hold_force,f'Condensation force = {np.round(hold_force,2)}',color = 'blue',rotation=00,size='medium')

save_fig = False
if save_fig:
    plt.savefig(filename + ".png", bbox_inches='tight')
    plt.close(fig)
