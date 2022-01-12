from bokeh.models import HoverTool, ColumnDataSource
from bokeh.plotting import figure
from bokeh.io import save, output_file
import sys
from operator import methodcaller
#The script is used to generate depth plot from samtools depth <bam> >> output.txt file 
#Can be used as stand-alone command-line tool or in a pipeline
#Integrated into process_fastq.sh

#Plot dimensions & command-line arguments
x_label = "Genome position"
y_label = 'Coverage depth'
plot_width = 1750
plot_height = 500
variable_list = sys.argv
path_to_depth_report_file, output_file_path = variable_list[1], variable_list[2]


def bar_plot(depth_dict):
    source = ColumnDataSource(data=depth_dict) #data provided as dict of lists
    p=figure(
        x_axis_label = x_label, #set label for x
        y_axis_label = y_label, #set label for y
        plot_width=plot_width, #set picture width in pixels
        plot_height=plot_height #set picture height in pixels
    )
    
    # #create hover tool with 2 tooltips, attached to the right of the bar
    hover = HoverTool(
    tooltips=[(x_label, "@{Genome position}"), (y_label, "@{Coverage depth}")],
    attachment = 'right'
    )
    
    p.xaxis.major_label_orientation = "vertical" #format label on x to orient vertically
    p.add_tools(hover) #add hover to the plot
    p.vbar(x='Genome position', top='Coverage depth',  width=0.5, bottom=0, fill_color=('green'), source=source, name = 'vbar')
    # #initiate bar plot object
    output_file(output_file_path)
    save(p, output_file_path)
    # show(p) #open file in browser

#Reading samtools depth report and converting it to {genome_position:coverage} dict
with open(path_to_depth_report_file, "r+") as depth_file:
    data_dict = dict(list(map(lambda x:list(map(int,x[1:])), list(map(methodcaller("split", "\t"), depth_file.read().split('\n')))))[:-1])
    depth_dict = {'Genome position':list(data_dict.keys()), 'Coverage depth':list(data_dict.values())}
    # print(list(depth_dict.items())[0:50])

#Generating html file with coverage report 
bar_plot(depth_dict)