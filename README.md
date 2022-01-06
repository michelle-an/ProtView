# ProtView

ProtView is an easy-to-use graphical sequence visualization, analysis, and data preparation tool that can compare .fasta formatted datasets using calculated biochemical features without relying on sequence homology based approaches. Graphs are made using Plotly.


I included a requirements.text file to easily manage dependancies
pip install -U pip
pip install -r requiremets.txt


Setup instructions to make ProtView easy to launch for mac OS X users:
1) Make sure you have python3 installed on your computer
	where python3   or   where python
2) Update the first line of the ProtView.command file after the #! to match the location of your python instalation.
3) Install python dependencies
	pip install -U pip
	pip install -r requirements.txt
4) Update user permissions
	chmod u+x ProtView.command
5) You can now double click on the ProtView.command file to run without needing to use Terminal.
