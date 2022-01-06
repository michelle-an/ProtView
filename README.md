# ProtView

ProtView is an easy-to-use graphical sequence visualization, analysis, and data preparation tool that can compare .fasta formatted datasets using calculated biochemical features without relying on sequence homology based approaches. Graphs are made using Plotly.
<br><hr><br>
I included a requirements.txt file to easily manage dependancies
<ul>
  <li>pip install -U pip</li>
  <li>pip install -r requiremets.txt</li>
</ul>
<br><hr><br>
Instructions to make ProtView easy to launch for mac OS X users:
<ol>
 	<li>Make sure you have python3 installed on your computer</li>
		<br>
		&emsp;
		<code>where python3</code>
		<br><br>
  	<li>Update the first line of the ProtView.command file after the #! to match the location of your python instalation</li>
		<br>
	<li>Install python dependencies</li>
		<br>
		&emsp;
		<code>pip install -U pip</code>
		<br>
		&emsp;
		<code>pip install -r requirements.txt</code>
		<br><br>
	<li>Update user permissions</li>
		<br>
		&emsp;
		<code>chmod u+x ProtView.command</code>
		<br><br>
	<li>You can now double click on the ProtView.command file to run without needing to use Terminal</li>
</ol>
