"""Pipeline configuration."""

from wordcloud import STOPWORDS

# Path names to files
_PROJ_HOME = "/gpfs/projects/bsc08/shared_projects/scGRN_analysis"
_DATA_HOME = f"{_PROJ_HOME}/Data_home/res/covid_19"
_META_FILE = f"{_PROJ_HOME}/Data_home/data/GSE145926_RAW/metadata.tsv"

# General purpose constants
# ggplot2 default colors
_GGPLOT_COLORS = dict(
    green="#39B600",
    yellow="#D89000",
    red="#F8766D",
    blue="#00B0F6",
    purple="#9590FF",
    cyan="#00BFC4",
    pink="E76BF3",
    light_pink="#FF62BC",
    saturated_green="#00BF7D",
)
# Setting up constant parameters
_SEED = 42  # random seed
_NODE_SIZE = 1200  # the node size when plotting
_ALPHA = 0.5  # transparency of the nodes/edges when plotting

_STOPWORDS = STOPWORDS.union(
    {
        "regulation",
        "activity",
        "positive",
        "negative",
        "catabolic",
        "process",
        "protein",
        "complex",
        "binding",
        "response",
        "gene",
        "genes",
        "encoding",
        "defining",
        "GeneID",
        "regulated",
    }
)
