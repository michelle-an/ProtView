# ProtView

ProtView is an easy-to-use graphical sequence visualization, analysis, and data preparation tool that can compare .fasta formatted datasets using calculated biochemical features without relying on sequence homology based approaches. Graphs are made using Plotly.
<br><hr>
### How to use ProtView

After launching ProtView a GUI will appear to guide the user through the setup process

<img width="700" alt="Screen Shot 2022-01-05 at 8 16 46 PM" src="https://user-images.githubusercontent.com/53797248/148327375-1d038e3d-6f17-492e-ad9e-d96ad3ad8158.png">

<ol>
	<li>
		Start by giving the project a name to easily find your results later.
	</li>
	<li>
		Next browse for .fasta formatted files using the Browse buttons on the left. You can optionally add a name for these files in the blank Class Name text field next to each button.
	</li>
	<li>
		Check or uncheck any of the boxes on the bottom to match which outputs you want.
	</li>
	<li>
		Select an output location to save your work using the Browse button at the bottom.
	</li>
	<li>
		Select Ok to begin, or press Cancel to exit.
	</li>
</ol>

<hr>

### Outputs

ProtView generates Plotly graphs to compare the distributions of 42 biochemical features between each of the class datasets that are uploaded.

A single joy plot is generated for each of the 42 feature

![joy plot](https://user-images.githubusercontent.com/53797248/148328033-f1b70a0b-ef55-40ab-a820-2b3cca74bc5a.png)

And statistical plots are generated to easily identify correlations in the data
	
![t-values](https://user-images.githubusercontent.com/53797248/148328202-de977189-ece7-4148-93c0-e15e8488a53c.png)

<br><hr>
### Instructions to make ProtView easy to launch for mac OS X users:
<ol>
 	<li>Make sure you have python3 installed on your computer</li>
		<br>
		&emsp;
		<code>
			where python3
		</code>
		<br><br>
  	<li>Update the first line of the ProtView.command file after the #! to match the location of your python instalation</li>
		<br>
	<li>Install python dependencies</li>
		<br>
		&emsp;
		<code>
			pip install -U pip
		</code>
		<br>
		&emsp;
		<code>
			pip install -r requirements.txt
		</code>
		<br><br>
	<li>Update user permissions</li>
		<br>
		&emsp;
		<code>
			chmod u+x ProtView.command
		</code>
		<br><br>
	<li>You can now double click on the ProtView.command file to run without needing to use Terminal</li>
</ol>
