import lumicks.pylake as lk
import matplotlib.pyplot as plt
import DNA
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

#Globally set the matplitlib font
plt.rcParams.update({'font.size': 18})


def plot_trace(trace):
    
    # Extract data from trace
    path = trace["Path"]
    filename = trace["Filename"]
    
    #filename = "20230317-130624 Kymograph 230317 Tether 1 Bleaching Kymo.h5"
    print(f"Analyzing file: {filename}")
    file = lk.File(f"{path}/{filename}")

    #Set up first figure with 3 subplots showing the blue and green kymographs
    # And the force-distance curve
    fig = plt.figure(1, layout="constrained")
    fig.clf()
    fig.suptitle(filename, fontsize = 16)
    #fig, axs = plt.subplots(4,1,constrained_layout=True, gridspec_kw={'height_ratios': [2, 2, 1, 1]})
    fig.subplots(4,1, gridspec_kw={'height_ratios': [2, 2, 1, 1]})
    
    #Make some adjustment to the spacing between subplots
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.0,
                        hspace=0.5)
    #Adjust the size of the subplots
    fig.set_figheight(15)
    fig.set_figwidth(15)
    
    #Plot the kymograph
    axes = plt.subplot(4,1,1)
    plt.title("")
    _, kymo = file.kymos.popitem()
    #Set the min and max photon count for the kymograph
    min_photon_count = 0
    max_photon_count = 50
    adjustment = lk.ColorAdjustment(min_photon_count, max_photon_count, mode="absolute")
    img = kymo.plot(channel="blue", adjustment = adjustment, aspect='auto', interpolation = 'none')
    
    #Aspect='auto'  this is a key parameter to stretch the kymo
    plt.xlabel("")
    #Add the image scalebar
    img_scale = fig.colorbar(img)
    img_scale.set_label("Photon Counter [Counts]")
    
    axes_img = plt.subplot(4,1,2)
    _, kymo = file.kymos.popitem()
    #Set the min and max photon count for the kymograph
    min_photon_count = 0
    max_photon_count = 50
    adjustment = lk.ColorAdjustment(min_photon_count, max_photon_count, mode="absolute")
    img = kymo.plot(channel="green", adjustment = adjustment, aspect='auto', interpolation = 'none')
    plt.title("")
    plt.xlabel("")
    #Add the image scalebar
    img_scale = fig.colorbar(img)
    img_scale.set_label("Photon Counter [Counts]")
    #Plot the force versus time
    axes = plt.subplot(4,1,3) #, sharex=axes
    axes.clear()
    #sharex=axes sends the axes from the above kymograph to this plot.  This
    # linkes the two plots x-axis together.
    line_timestamp_ranges = kymo.line_timestamp_ranges()
    force = file.force2x
    downsampled_force = force.downsampled_over(line_timestamp_ranges)
    #force.plot(label="high frequency") #Don't plot the HF force
    #Do plot the downsamples force.
    downsampled_force.plot(label="downsampled to kymo rate")
    #This line finds the duration of the kymographs.  It is used to set the 
    # x-axis limits below.
    time_end  = (line_timestamp_ranges[-1][0] - line_timestamp_ranges[0][0])*(10**-9)
    axes.sharex(axes_img)
    axes.set_ylim(0,2)
    #Plot the force versus extension
    axes = plt.subplot(4,1,4)
    #Plot the thoery first so it is underneath the data
    Lp = 50
    K0 = 1200
    kbT = 4.09
    bp_lambda = 48502
    t_force = np.arange(0,60,0.02)
    t_ext = DNA.extension_MMS(t_force,Lp,K0,kbT)*bp_lambda*0.000334
    plt.plot(t_ext,t_force, color='grey', linestyle='dashed')
    #Plot the force extension data
    extension = file.distance1
    downsampled_extension = extension.downsampled_over(line_timestamp_ranges)
    #I had a wierd issue where the downsampled force and downsampled distance 
    # had a length different by 1?  I added this code to crop the longer array.
    length = min(len(downsampled_extension.data),len(downsampled_force.data))-1
    plt.plot(downsampled_extension.data[:length],downsampled_force.data[:length])
    plt.ylabel("Force (pN)")
    plt.xlabel("Extension (um)")
    #Hardcoded the axis limits.  They can be changed as needed.
    axes.set_xlim(0,18)
    axes.set_ylim(0,2)
    

    #Save the figure to the current path with the same filename + .png
    #plt.savefig(f"{path}/{filename}.png")
    #Show the figure
    plt.savefig(filename + "_kymo.png")
    plt.show()
    
    ## End of figure 1
    

    #######################################################################################
    #######################################################################################
    ####################################################################################### 
    # Plot the kymo for figure
    ## Start figure 2
    blue_img = kymo.get_image("blue")
    green_img = kymo.get_image("green")
    line_time_seconds = kymo.line_time_seconds #time in of a line in seconds.
    pixel_size = kymo.pixelsize_um[0]
    
    try:
        time_start = float( trace["Start"])
    except KeyError:
        time_start = 0
    
    try:
        time_end = float(trace["End"])
    except KeyError:
        time_end = np.size(blue_img[0,:])*line_time_seconds
    frame_start = round(time_start/line_time_seconds)
    frame_end = round(time_end/line_time_seconds)
    
    # Blue and green excitation is interlaced in kymograph to reduce color channel crosstalk.
    # This code de-interlaces the kymograph.
    blue_img_crop = blue_img[:,frame_start:frame_end]
    green_img_crop = green_img[:,frame_start:frame_end]
    frame_select = np.mean(blue_img_crop[0:10,:],axis=0) > np.mean(blue_img_crop[0:10,:]) # Use the blue bead to select frames.
    b_frame_select = frame_select
    g_frame_select = [not item for item in frame_select]

    index_start = 0
    index_end = int(round(10/pixel_size))
    
    blue_img_crop = blue_img[index_start:index_end,frame_start:frame_end]
    green_img_crop = green_img[index_start:index_end,frame_start:frame_end]
    
    frame_select = np.mean(blue_img_crop[0:10,:],axis=0) > np.mean(blue_img_crop[0:10,:]) # Use the blue bead to select frames.
    
    
    b_frame_select = frame_select
    g_frame_select = [not item for item in frame_select]
    #g_frame_select = frame_select
    
    #b_frame_select = slice(None)
    #g_frame_select = slice(None)
    (numrows,numcols) = np.shape(blue_img_crop)
    
    extent = (0, numcols*line_time_seconds, numrows*pixel_size, 0)
    
    fig = plt.figure(3,layout="constrained")
    fig.clf()
    fig.set_figheight(13.75)
    fig.set_figwidth(11.5)
    fig.suptitle(filename, fontsize = 16)

    try:
        b_max = float(trace["b_max"])
        g_max = float(trace["g_max"])
    except:
        print("Warning max intensity not found")
        b_max = 75 # max value for blue channel
        g_max = 75 # max value for green channel
    
    ax1 = plt.subplot(4,2,1)
    plt.title("Initial Kymo")
    # Plot the blue channel as red since it has better contrast
    b_cmap = LinearSegmentedColormap.from_list("b_cmap", ["black","red"])
    b_plt = plt.imshow(blue_img_crop[:,b_frame_select], cmap = b_cmap, vmin=0, vmax = b_max, aspect = 'auto', interpolation= 'none', extent=extent)
    plt.ylabel("Position (um)")
    plt.xlabel("Time (s)")
    
    ax2 = plt.subplot(4,2,2)
    plt.title("Initial Kymo")
    g_cmap = LinearSegmentedColormap.from_list("g_cmap", [[0,0,0],[0,1,0]])
    g_plt = plt.imshow(green_img_crop[:, g_frame_select], cmap = g_cmap, vmin=0, vmax = g_max, aspect = 'auto', interpolation= 'none', extent=extent)

    plt.xlabel("Time (s)")
    plt.ylabel("Position (um)")
    
    #Get image from end
    time_start = 90
    time_end = 110
    frame_start = round(time_start/line_time_seconds)
    frame_end = round(time_end/line_time_seconds)
    blue_img_crop = blue_img[index_start:index_end,frame_start:frame_end]
    green_img_crop = green_img[index_start:index_end,frame_start:frame_end]
    
    frame_select = np.mean(blue_img_crop[0:10,:],axis=0) > np.mean(blue_img_crop[0:10,:]) # Use the blue bead to select frames.
    b_frame_select = frame_select
    g_frame_select = [not item for item in frame_select]
    
    
    ax3 = plt.subplot(4,2,3)
    plt.title("Post Kymo")
    # Plot the blue channel as red since it has better contrast
    b_cmap = LinearSegmentedColormap.from_list("b_cmap", ["black","red"])
    plt.imshow(blue_img_crop[:,b_frame_select], cmap = b_cmap, vmin=0, vmax = b_max, aspect = 'auto', interpolation= 'none', extent=extent)
    plt.xlabel("Time (s)")
    plt.ylabel("Position (um)")
    
    ax4 = plt.subplot(4,2,4)
    plt.title("Post Kymo")
    g_cmap = LinearSegmentedColormap.from_list("g_cmap", [[0,0,0],[0,1,0]])
    plt.imshow(green_img_crop[:, g_frame_select], cmap = g_cmap, vmin=0, vmax = g_max, aspect = 'auto', interpolation= 'none', extent=extent)
    plt.xlabel("Time (s)")
    plt.ylabel("Position (um)")
    
    # Get green and blue image scan
    ax5 = plt.subplot(4,2,5)
    scans = list(file.scans)
    
    try:
        img_left = float(trace["img_left"])
    except:
        img_left = 2
    img_top = 2.4
    img_height = 7
    img_width = 11
    
    scan = file.scans[scans[0]]
    b_img = scan.get_image("blue")
    pixel_size = scan.pixelsize_um[0]
    (numrows,numcols) = np.shape(b_img)
    b_img_crop = b_img[int(img_top/pixel_size):int((img_height+img_top)/pixel_size),int(img_left/pixel_size):int((img_width+img_left)/pixel_size)]
    (numrows,numcols) = np.shape(b_img_crop)
    extent = (0, numcols*pixel_size, numrows*pixel_size, 0)
    plt.imshow(b_img_crop, cmap = b_cmap, vmin=0, vmax = b_max, aspect = "equal", interpolation= 'none', extent=extent)
    plt.ylabel("Y Position (um)")
    plt.xlabel("X Position (um)")
   
    ax6 = plt.subplot(4,2,6)
    scan = file.scans[scans[1]]
    g_img = scan.get_image("green")
    pixel_size = scan.pixelsize_um[0]
    (numrows,numcols) = np.shape(g_img)
    
    g_img_crop = g_img[int(img_top/pixel_size):int((img_height+img_top)/pixel_size),int(img_left/pixel_size):int((img_width+img_left)/pixel_size)]
    (numrows,numcols) = np.shape(g_img_crop)
    extent = (0, numcols*pixel_size, numrows*pixel_size, 0)
    plt.imshow(g_img_crop, cmap = g_cmap, vmin=0, vmax = g_max, aspect = "equal", interpolation= 'none', extent=extent) #
    plt.ylabel("Y Position (um)")
    plt.xlabel("X Position (um)")
    
    # Merge two colors
    ax7 = plt.subplot(4,2,7)
    r_img = np.zeros(g_img_crop.shape)
    color_cmap = LinearSegmentedColormap.from_list("color_cmap", [[0,0,0],[1,1,0]])
    plt.imshow(np.dstack((r_img,g_img_crop/g_max,b_img_crop/b_max)), cmap = color_cmap, aspect = "equal", interpolation= 'none', extent=extent) 

    img_scale = plt.colorbar(b_plt, ax= [ax1,ax3], shrink = 0.5)
    img_scale.set_label("Blue Photon Counter [Counts]")
    
    img_scale = plt.colorbar(g_plt, ax= [ax2,ax4], shrink = 0.5)
    img_scale.set_label("Green Photon Counter [Counts]")

    plt.savefig(filename+"_example.png")
    plt.show()
    ## End of plot_trace fucntion    



###############################################################################
# Script to analyze a list of traces
###############################################################################

def make_dict(keys, items):
    return {key:item for (key,item) in zip(keys,items)}

# Load a trace description file
with open("traces_for 2scan_with max intenisty values.txt") as f:
    keys = [item.strip() for item in f.readline().split("\t")]
    traces = [make_dict(keys, [item.strip() for item in line.split("\t")]) for line in f.readlines()]

#Plot a single trace
# print(traces)
traces = [traces[3]]
# Loop over traces, make and save image figures.
[plot_trace(trace) for trace in traces]
