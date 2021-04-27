# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:27:18 2021

@author: Trevor.Grout
"""

from pathlib import Path
import pandas as pd
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


#########################################################################
#Create scatter plot
#########################################################################
def scatterplot(dataframe, x_field, y_field, hue_field, title_text, stats_text=False, annotate = False, dest_file = False):
    '''
    Create boxplots.

    Parameters
    ----------
    dataframe : DataFrame
        Pandas dataframe data to be plotted.
    x_field : STR
        Field to use for x-axis (Assumes FIM 2)
    y_field : STR
        Field to use for the y-axis (Assumes FIM 3)
    title_text : STR
        Text for plot title.
    stats_text : STR or BOOL
        Text for stats to place on chart. Default is false (no stats printed)
    dest_file : STR or BOOL, optional
        If STR provide the full path to the figure to be saved. If False
        no plot is saved to disk. The default is False.

    Returnsy
    -------
    fig : MATPLOTLIB
        Plot.

    '''

    #initialize plot
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 12))

    #Use seaborn to plot the boxplot
    axes=sns.scatterplot(data=dataframe, x=x_field, y=y_field, hue = hue_field, color = 'black', s = 150)
    #sns.regplot(x=x_field, y=y_field, data=dataframe, scatter=False, ax=axes)
    #Set xticks and yticks and background horizontal line.
    ymax = dataframe[y_field].max()
    xmax = dataframe[x_field].max()
    axes.set(ylim=(0.0,ymax),yticks = np.arange(0,ymax,round(ymax/10,1)))
    axes.set(xlim=(0.0,xmax),xticks = np.arange(0,xmax,round(xmax/10,1)))
    axes.grid(b=True, which='major', axis='both')

    #Set sizes of ticks and legend.
    axes.tick_params(labelsize = 'xx-large')

    #Define y axis label and x axis label.
    axes.set_ylabel(f'{y_field.replace("_"," ")}',fontsize='xx-large',weight = 'bold')
    axes.set_xlabel(f'{x_field.replace("_"," ")}',fontsize='xx-large',weight = 'bold')

    #set title of plot
    axes.set_title(f'{title_text}',fontsize=20, weight = 'bold')
    
    #rename legend labels to the simplified labels.
    axes.legend(markerscale = 2.5, fontsize = 25)


    if annotate:
        #Set text for labels
        box_props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        textbox_str = 'Target Better'
        axes.text(0.3, 0.6, textbox_str, transform=axes.transAxes, fontsize=32, color = 'gray', fontweight = 'bold', verticalalignment='top', bbox=box_props, rotation = 35, rotation_mode = 'anchor')
        textbox_str = 'Baseline Better'
        axes.text(0.5, 0.2, textbox_str, transform=axes.transAxes, fontsize=32, color = 'gray', fontweight = 'bold', verticalalignment='top', bbox=box_props, rotation = 35, rotation_mode = 'anchor')

    if stats_text:
        #Add statistics textbox
        axes.text(0.01, 0.80, stats_text, transform=axes.transAxes, fontsize=24, verticalalignment='top', bbox=box_props)

    #If figure to be saved to disk, then do so, otherwise return fig
    if dest_file:
        fig.savefig(dest_file)
        plt.close(fig)
    else:
        return fig
    

    
eval_plots_file = '/path/to/fim_performance_points.shp'
usgs_gages = '/path/to/usgs_gages.gpkg'
sierra_file = 'path/to/sierra/test/output/agg_nwm_recurr_flow_elev_stats_location_id_recurr_interval.csv'



fim_gdf = gpd.read_file(eval_plots_file)
usgs_gdf = gpd.read_file(usgs_gages, dtype={'location_id':str})
sierra_df = pd.read_csv(sierra_file, dtype={'location_id':str})

usgs_gdf['nws_lid'] = usgs_gdf['nws_lid'].str.lower()


#Create spatial layer of Sierra mean_abs_y_diff_ft across all events for a site.
agg = sierra_df.groupby('location_id')['mean_abs_y_diff_ft'].agg(['mean','size']).reset_index()
gages = usgs_gdf.merge(agg[['location_id', 'mean']], on = 'location_id').dropna(subset=['mean'])
gages.to_file(Path(eval_plots_file).parent / 'joined_all_sites.shp')







all_gdf=gpd.GeoDataFrame()

for category in ['action','minor','moderate','major']:
    cat_fim = fim_gdf.query(f'magnitude == "{category}"').drop(columns=['geometry'])
    cat_sierra = sierra_df.query(f'recurr_interval == "{category}"')
    
    joined = usgs_gdf.merge(cat_fim[['nws_lid','version','magnitude','source','CSI']], on = 'nws_lid')
    joined = joined.merge(cat_sierra[['location_id','mean_abs_y_diff_ft']], on = 'location_id')

    #joined = joined.query('source == "usgs"')
    #joined = joined.query(f'mean_abs_y_diff_ft < {joined.mean_abs_y_diff_ft.quantile(0.85)}')  
    joined.rename(columns={'mean_abs_y_diff_ft':'mean_absolute_difference_(ft)'}, inplace = True)

    today = datetime.today().strftime('%Y-%m-%d')
    version = joined.version.unique()[0].replace('fim_','')
        
    dest_file = Path(eval_plots_file).parent / f'{category}.png'
    plot = scatterplot(joined,'mean_absolute_difference_(ft)', 'CSI', 'source', f'CSI vs Mean Abs Diff (ft)\n{category.title()} Category, FIM Version {version}\nGID, NWC, OWP\n{today}', False, False, dest_file)
    
    all_gdf = all_gdf.append(joined)
