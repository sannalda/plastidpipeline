import os
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.express as px


def setup_logger(log_file_path):
    loggerPP = logging.getLogger("PlastidPipeline_%s" %snakemake.wildcards["sample"])
    loggerPP.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
    console_handler.setFormatter(formatter)
    loggerPP.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    loggerPP.addHandler(file_handler)

    return loggerPP

#loggerPP = setup_logger(os.path.join(snakemake.config["workdir"],snakemake.log[0]))


error_log_file = os.path.join(snakemake.config["workdir"],snakemake.log[0])
sys.stderr = sys.stdout = open(error_log_file, "w+")

# Adapted from StackOverflow user vestland (thanks!)
def highLights(fig, variable, level, mode, fillcolor, layer):
    """
    Set a specified color as background for given
    levels of a specified variable using a shape.
    
    Keyword arguments:
    ==================
    fig -- plotly figure
    variable -- column name in a pandas dataframe
    level -- int or float
    mode -- set threshold above or below
    fillcolor -- any color type that plotly can handle
    layer -- position of shape in plotly figure, like "below"
    """
    
    if mode == 'above':
        m = df[variable].gt(level)
    
    if mode == 'below':
        m = df[variable].lt(level)
        
    df1 = df[m].groupby((~m).cumsum())['pos'].agg(['first','last'])

    for index, row in df1.iterrows():
        fig.add_shape(type="rect",
                        xref="x",
                        yref="paper",
                        x0=row['first'],
                        y0=0,
                        x1=row['last'],
                        y1=1,
                        line=dict(color="rgba(0,0,0,0)",width=3,),
                        fillcolor=fillcolor,
                        layer=layer) 
        print("Range of loci with number of reads %s %d: %d-%d" %(mode,level,row['first'],row['last']))
        #loggerPP.debug("Range of loci with number of reads %s %d: %d-%d" %(mode,level,row['first'],row['last']))
    return(fig)

df = pd.read_csv(os.path.join(snakemake.config["workdir"],snakemake.input[0]),sep = "\t",names=["sample","pos","reads"])
fig = px.line(df, x="pos", y="reads")
fig.add_hline(y=snakemake.config["Backmapping"]["MinReadCoverageWarning"])
fig = highLights(fig = fig, variable = 'reads', level = snakemake.config["Backmapping"]["MinReadCoverageWarning"], mode = 'below',
               fillcolor = 'rgba(255,0,0,0.75)', layer = 'below')
fig.update_layout(title='Number of Reads Backmapped to Assembled Genome, Sample %s' %snakemake.wildcards["sample"],
                   xaxis_title='Position',
                   yaxis_title='Number of Reads')
#fig.show()
fig.write_html(os.path.join(snakemake.config["workdir"],snakemake.output[0]))