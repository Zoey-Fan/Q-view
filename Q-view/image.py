"""
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  Curve line-of-sight visibility overlay map.
 """

###merge all csvs 
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

csv_dir='/home/zoe/Documents/Experiment-viewshed/Q-view/Q-view/csv-curve'
csv_files=glob.glob(os.path.join(csv_dir,'*csv'))
output_csv='merged_data.csv'
merged_data=pd.DataFrame(0,index=range(500),columns=range(500))
# Screen coordinates and tile size
Min_X=0
Min_Y=0
tile_size=500

# Loop through all CSV files found
for file in csv_files:
    filename_without_ext = os.path.splitext(file)[0] 
    # Split the filename using underscores to get coordinates 
    parts = filename_without_ext.split('_')
    first_number = parts[1]  
    second_number = parts[2]  
    X = int(first_number)  
    Y = int(second_number)  
    # Calculate the top-left corner to place the data (with a margin of 'r')
    radius=100
    interestRange_x=max(Min_X,X-radius)
    if interestRange_x>0:
        interestRange_x=interestRange_x-1
    interestRange_y=max(Min_Y,Y-radius)
    if interestRange_y>0:
        interestRange_y=interestRange_y-1

    df=pd.read_csv(file,header=None)
    # Flip the DataFrame vertically (reverse the rows)
    df_flipped=df.iloc[::-1]

    start_row,start_col=interestRange_y,interestRange_x
    if start_row+df.shape[0]>tile_size or start_col+df.shape[1]>tile_size:
        raise ValueError("The small file is too large to fit into the large file.")  
    # Add the data from the small file to the corresponding position in the large file
    merged_data.iloc[start_row:start_row+df.shape[0],start_col:start_col+df.shape[1]]+=df.values

merged_data.to_csv('merged_file.csv', header=False, index=False)

###printing into png
#0---write,1---green
assert merged_data.shape==(500,500),"Data is not 500*500"
fig,ax=plt.subplots()

#Define a custom color map for the heatmap
color = ['#C1FFC1', '#B4EEB4', '#4EEE94', '#00CD66','#008B00'] 
colors = LinearSegmentedColormap.from_list("GreenGradient", color)
merged_data=merged_data.replace(0,np.nan)

im=ax.imshow(merged_data,cmap=colors,interpolation='nearest',vmin=1,vmax=merged_data.max().max(),origin='lower')
# Set x and y ticks for the heatmap
xticks=np.arange(0,501,50)
yticks=np.arange(0,501,50)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
plt.xlabel('pixel')
plt.ylabel('pixel')
# Add a color bar to the heatmap
cbar=fig.colorbar(im,shrink=0.6,location='right',pad=0.05)
cbar.ax.tick_params(labelsize=8)
# Save the heatmap as a PNG file
plt.savefig('curve_image',dpi=300,bbox_inches='tight',pad_inches=0)
plt.show()


