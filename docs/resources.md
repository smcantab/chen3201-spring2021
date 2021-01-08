# Hello world

## Python

There are many great resources for learning Python. Here are some of my and my students' favourites.

### Conda/Anaconda/pip
For a local installation of Python with many of the Data Science libraries you may want to use, I recommend installing Conda/Anaconda.

- [Anaconda](https://www.anaconda.com/products/individual)
    - [Installers](https://www.anaconda.com/products/individual)
    - [Documentation](https://docs.anaconda.com/anaconda/)
    - [Available ML packages](https://www.anaconda.com/open-source)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) is a package installer for Anaconda.
- [pip](https://pip.pypa.io/en/stable/quickstart/) is the Python package installer, it works similarly to conda.

## Integrated Development Environment (IDE)
There are many great IDEs to write Python code, I list below some of the most popular and easier to use.

- [Spyder](https://www.spyder-ide.org/), ships with [Anaconda][1]. It was developed specifically for scientific computing and it has a Matlab-like feel. Unless you already use another IDE, or have a strong preference, go with this.
- [PyCharm](https://www.jetbrains.com/pycharm/). You can get the professional edition for free registering with your UMN account.
- [Atom](https://atom.io/) is what I use. It is highly customizable, but needs some work to be set up as a Python IDE, _e.g._ see this [blog](https://medium.com/issuehunt/20-atom-plug-ins-for-python-development-d6b10f8fa33e).
- [Sublime Text](https://www.sublimetext.com/) is similar to Atom, to set it up as a Python IDE see for instance this [blog](https://realpython.com/setting-up-sublime-text-3-for-full-stack-python-development/).

### Guides/Tutorials

- [Matlab vs. Python](https://realpython.com/matlab-vs-python/)
- [Real Python Tutorials](https://realpython.com/)
- [Python like you mean it](https://www.pythonlikeyoumeanit.com/index.html)
- [Python for beginners](https://www.python.org/about/gettingstarted/)
- [Official Python 3 Tutorial](https://docs.python.org/3/tutorial/)
-  [Geeks for Geeks getting started with Python](https://www.geeksforgeeks.org/python-programming-language/?ref=lbp)
- [Problem Solving with Algorithms and Data Structures using Python](https://runestone.academy/runestone/books/published/pythonds/index.html). This is a great online book. It teaches data structures and algorithms while also teaching Python in quite some depth. It’s really easy to read and it will teach you all the basics to become a good (Python) programmer.
- [Moving to Python from other languages](https://wiki.python.org/moin/MovingToPythonFromOtherLanguages)

### IPython/Jupyter notebooks
- [Jupyter](https://jupyter.org/)/[Jupyter Lab](https://jupyterlab.readthedocs.io/en/stable/) are shipped with [Anaconda][1]
- [nbviewer](https://nbviewer.jupyter.org/) is another way of sharing Jupyter notebooks
- [Google Colaboratory (Colab)](https://colab.research.google.com/notebooks/intro.ipynb) enables collaborative editing of ipython notebooks and has its own computational resource (even GPUs! [How to enable GPUs and other Colab tutorials](https://www.tutorialspoint.com/google_colab/google_colab_using_free_gpu.htm))
- [Binder](https://mybinder.org/) is a tool to turn your GitHub repository in a "binder" of iPython notebooks. It is hosted on JupyterHub, but no easy/free GPUs I'm afraid!
- [MSI Jupyter Notebooks service](https://www.msi.umn.edu/support/faq/how-do-i-get-started-jupyter-notebooks). If you need more computational resources I suggest you use the MSI Jupyter Notebook service. I am negotiating some compute time for you to use on MSI as part of this course. Alternatively, you can register as part of a PI's research group, if they are willing to lend you some of their MSI compute time.

## Libraries
This is a nonexhaustive list of some of the better known Python libraries in this area, or less known but that I know of and like. Many of these libraries are shipped with the [Anaconda][1] installation.

[1]:resources.md#condaanaconda

### Scientific Computing
- [NumPy](https://numpy.org/devdocs/user/quickstart.html)
    - [Geeks for Geeks Numpy in Python](https://www.geeksforgeeks.org/numpy-in-python-set-1-introduction/)
- [SciPy](https://www.scipy.org/)
    - [Geeks for Geeks Data Analysis with Scipy](https://www.geeksforgeeks.org/data-analysis-with-scipy/)

### ML & Data Science
- [Pandas](https://pandas.pydata.org/docs/index.html)
- [scikit-learn](https://scikit-learn.org/stable/)
    - [Choosing the right scikit-learn estimator](https://scikit-learn.org/stable/tutorial/machine_learning_map/index.html)
- [PyTorch](https://pytorch.org/)
    - [Tutorials](https://pytorch.org/tutorials/beginner/deep_learning_60min_blitz.html) (60 minutes blitz)
- [Tensorflow](https://www.tensorflow.org/)
    - [Tutorials](https://www.tensorflow.org/tutorials)
    - [Keras](https://keras.io/getting_started/intro_to_keras_for_researchers/)
    - [Sonnet](https://github.com/deepmind/sonnet)
- [MXNet](https://mxnet.apache.org/versions/1.6/)
    - [Tutorials](https://mxnet.apache.org/versions/1.6/api/python/docs/tutorials/)
    - [Gluon](https://mxnet.apache.org/versions/1.6/api/python/docs/tutorials/getting-started/crash-course/index.html) (60 minutes blitz)
- [JAX](https://github.com/google/jax)
    - [Quickstart](https://jax.readthedocs.io/en/latest/notebooks/quickstart.html)
- [ONNX](https://github.com/onnx/onnx)
- [PyTorch Geometric](https://github.com/rusty1s/pytorch_geometric) (there are more)
- [Deep Graph Library (DGL)](https://github.com/dmlc/dgl)

### ML for Chemistry and Materials

- [DeepChem](https://deepchem.io/)
    - [Book: _Deep learning for the life sciences: applying deep learning to genomics, microscopy, drug discovery, and more_](https://www.oreilly.com/library/view/deep-learning-for/9781492039822/)
    -  [Tutorials](https://github.com/deepchem/deepchem/tree/master/examples/tutorials)
- [JAX M.D.](https://github.com/google/jax-md)
    - [Paper](https://arxiv.org/abs/1912.04232)
    - [Quickstart](https://github.com/google/jax-md#getting-started)
- [ML4Chem](https://github.com/muammar/ml4chem)
- [DScribe](https://singroup.github.io/dscribe/latest/)
- [KLIFF](https://kliff.readthedocs.io/en/stable/)
- [TorchANI](https://github.com/aiqm/torchani)
    - [Paper](https://pubs.acs.org/doi/10.1021/acs.jcim.0c00451)
- [SchnetPack](https://github.com/atomistic-machine-learning/schnetpack)
    - [Paper](https://pubs.acs.org/doi/10.1021/acs.jctc.8b00908)
- [DeepPMD](https://github.com/deepmodeling/deepmd-kit)
    - [Manual](http://www.deepmd.org/deepmd-kit-manual/)
- [ASAP](https://github.com/BingqingCheng/ASAP)
    - [Paper](https://pubs.acs.org/doi/full/10.1021/acs.accounts.0c00403)
- [AMP](https://amp.readthedocs.io/en/latest/)
    - [Paper](https://www.sciencedirect.com/science/article/pii/S0010465516301266?via%3Dihub)

### Visualization
- [Seaborn](https://seaborn.pydata.org/)
    - [Example gallery](https://seaborn.pydata.org/examples/index.html)
- [Matplotlib](https://matplotlib.org/)
    - [Geeks for Geeks Introduction to Matplotlib](https://www.geeksforgeeks.org/python-introduction-matplotlib/)
- [Plotly Dash](https://plotly.com/dash/)
    - [Dash gallery](https://dash-gallery.plotly.host/Portal/)
    - [GitHub page](https://github.com/plotly/plotly.py)

### Cheat Sheets
- [Python for Data Science Cheat Sheet](https://s3.amazonaws.com/assets.datacamp.com/blog_assets/Numpy_Python_Cheat_Sheet.pdf)
- [Pandas Cheat Sheet](https://pandas.pydata.org/Pandas_Cheat_Sheet.pdf)
- [Matplotlib Cheat Sheet](https://github.com/matplotlib/cheatsheets)
- [Conda Cheat Sheet](https://conda.io/projects/conda/en/latest/_downloads/843d9e0198f2a193a3484886fa28163c/conda-cheatsheet.pdf)

## Data Sets
- [QM# and other data sets](http://quantum-machine.org/)
- [MoleculeNet](http://moleculenet.ai/datasets-1)
- [AFLOW](http://www.aflowlib.org/)
- [CMR](https://cmr.fysik.dtu.dk/)
- [OQMD](http://oqmd.org/)
- [QCArchive](https://qcarchive.molssi.org/)
- [kaggle _Predicting Molecular Properties_](https://www.kaggle.com/c/champs-scalar-coupling/overview)

## Standards
- [OPTIMADE](https://www.optimade.org/)
- [MolSSI](https://molssi.org/software-projects/)

## Git/GitHub
Git is a version control system that enable the collaborative development/storage of Software. GitHub is the most popular code hosting platform for version control and collaboration.

- [Git Installation](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [GitHub guides](https://guides.github.com/)
    - [Getting started](https://guides.github.com/activities/hello-world/)
    - [YouTube channel](https://www.youtube.com/githubguides)

## Stack Overflow
[Stack Overflow](https://stackoverflow.com/) is a great forum to find answers to your programming questions and ask new ones if nobody has aver asked the same question. Parsing Stack Overflow is one of the best ways of learning. If you have a programming question, before writing a post on [Piazza](https://piazza.com/umn/fall2020/chen5595) or contacting the [TA](https://www.cems.umn.edu/people/grads/rishabh-gupta), I suggest you check Stack Overflow.

## Other UMN computing Resources
- [UMN Virtual Online Linux Environment](https://vole.cse.umn.edu/) is a full CSE Labs Linux desktop delivered virtually through an application called FastX. This means you can do homework, programming and even use engineering software from your laptop. As you can connect to the internet, you can connect to a VOLE desktop.

## Tools for teams
- [Slack](https://slack.com/) is a communication platform for teams that integrates well with a large number of tools (_e.g._ Google Drive/Calendar, Zoom, GitHub, Workstreams.ai, Slite, etc.)
- [Workstreams.ai](https://workstreams.ai/) is a Kanban board-like service for task management
- [Slite](https://slite.com/) is an internal knowledge-base platform for teams. It has a Slack-like interface and it is great way or organizing/sharing notes/knowledge.

## Maths/Stats
In addition to the textbooks listed under [Organization/Books][1], you may find useful the following resources:

- [The Matrix Cookbook](matrixcookbook.pdf) is a quick reference guide for matrix identities, relations, etc.
- [d2l.ai Maths for Deep Learning](https://d2l.ai/chapter_appendix-mathematics-for-deep-learning/index.html) is an informal overview of some of the mathematics required for this course, most of which you know already or that we will cover in the first part of the course.
- [A Statistical Learning/Pattern Recognition Glossary](http://alumni.media.mit.edu/~tpminka/statlearn/glossary/) by Thomas Minka
- Datasets with Varied Appearance and Identical Statistics: **The datasaurus dataset**
    - [Paper](SameStats-DifferentGraphs.pdf)
    - [Blog](https://www.autodeskresearch.com/publications/samestats)

[1]: organization.md#books
