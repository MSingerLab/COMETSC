import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

project = 'COMET'
highlight_language = 'python'
copyright = '2019, Meromit Singer'
author = 'Conor Delaney'
# version and release are intended to provide a separation between 'X.Y' short
# version number ('version') and 'X.Y.Z-A.B' full release number
# ('release'). We don't need this separation, so set both to same value to
# prevent confusion.
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_favicon = '_static/COMET_favicon.ico'
html_short_title = 'COMET'
