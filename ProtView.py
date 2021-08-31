#!/Users/michellean/opt/anaconda3/bin/python

import os
import sys
import time
import datetime

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import progressbar
import PySimpleGUI as sg
import statsmodels.api
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from plotly.subplots import make_subplots
from sklearn.ensemble import RandomForestRegressor


def fileQC(files):
    if len(files) == 0:
        sys.exit(print("No input files provided"))

    bad_files = []

    for file in files:
        if not file.endswith(".fasta"):
            bad_files.append(file)
        else:
            with open(file, "r") as f:
                f = f.read()
            sequences_with_headers = [x.strip("\n") for x in f.split(">")][1:]
            sequences = [
                "\n".join(x.split("\n")[1:]).upper() for x in sequences_with_headers
            ]

            seq_num = 0
            for seq in sequences:
                seq_num += 1
                aas = set(seq)

                possible_values = [
                    "A",
                    "B",
                    "C",
                    "D",
                    "E",
                    "F",
                    "G",
                    "H",
                    "I",
                    "J",
                    "K",
                    "L",
                    "M",
                    "N",
                    "O",
                    "P",
                    "Q",
                    "R",
                    "S",
                    "T",
                    "U",
                    "V",
                    "W",
                    "Y",
                    "Z",
                    "X",
                    "*",
                    "-",
                    " ",
                    "\n",
                    "\t",
                ]

                for aa in aas:
                    if aa not in possible_values:
                        bad_files.append(f"{file}\t\tsequence #: {seq_num}")
                        break

    return bad_files


def simpleGUI():
    sg.theme("DarkBlue3")
    layout = [
        [
            sg.Text("Instructions:", size=(10, 4), font=("Helvetica", 16)),
            sg.Text(
                "1) Write project name.\n"
                "2) Select FASTA formatted files to evaluate.\n"
                "3) By default, glimpse will make all types of plots. Check the specific boxes if you want to skip any.\n"
                "4) Designate output folder.",
                size=(90, 5),
                font=("Helvetica", 16),
            ),
            # sg.Text("2)Select FASTA formatted files to evaluate.", size=(90, 1), font=("Helvetica", 16), ),
            # sg.Text("3)designate output folder.", size=(90, 1), font=("Helvetica", 16), )
        ],
        [
            sg.Text("Project Name", font=("Helvetica", 12)),
            sg.InputText(
                key="project_name",
                size=(50, 1),
                default_text=str(int(round(time.time(), 0))),
                font=("Helvetica", 12),
            ),
        ],
        [sg.T("")],
        [
            sg.Text(
                "Select files and input classification names (optional)",
                size=(100, 1),
                font=("Helvetica", 16),
            )
        ],
        [sg.T("")],
        [
            sg.FileBrowse(size=(10, 1), key="file1"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name1", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file2"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name2", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file3"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name3", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file4"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name4", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file5"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name5", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file6"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name6", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file7"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name7", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file8"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name8", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file9"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name9", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [
            sg.FileBrowse(size=(10, 1), key="file10"),
            sg.Text("Class Name: ", font=("Helvetica", 12)),
            sg.InputText(key="name10", size=(30, 1), font=("Helvetica", 12)),
            sg.Text("", size=(90, 1)),
        ],
        [sg.Text("")],
        [
            sg.Checkbox(
                "Show Graphs", default=False, key="show", font=("Helvetica", 12)
            ),
            sg.Checkbox(
                "Skip T-Values", default=False, key="skip_t", font=("Helvetica", 12)
            ),
            sg.Checkbox(
                "Skip RF", default=False, key="skip_rf", font=("Helvetica", 12)
            ),
            sg.Checkbox(
                "Skip Correlations", default=False, key="skip_c", font=("Helvetica", 12)
            ),
            sg.Checkbox(
                "Skip Joy Plots", default=False, key="skip_j", font=("Helvetica", 12)
            ),
            sg.Checkbox(
                "Save .csv of features",
                default=False,
                key="save",
                font=("Helvetica", 12),
            ),
        ],
        [sg.Text("Output Folder:", font=("Helvetica", 12))],
        [
            sg.FolderBrowse(size=(10, 1), key="out_browser"),
            sg.InputText(
                key="output",
                size=(50, 1),
                default_text=f"{os.getcwd()}",
                font=("Helvetica", 12),
            ),
        ],
        [
            sg.Button("Cancel", font=("Helvetica", 16)),
            sg.Button("Ok", font=("Helvetica", 16)),
        ],
    ]

    window = sg.Window("Glimpse", layout)

    while True:
        event, values = window.read()
        clean_values = {k: v for k, v in values.items() if v is not None}
        if event == sg.WIN_CLOSED or event == "Cancel":
            sys.exit(print("Stopped"))
        if event == "Ok":
            files = [v for k, v in values.items() if k.startswith("file")]
            names = [v for k, v in values.items() if k.startswith("name")]
            show = clean_values["show"]
            skip_t = clean_values["skip_t"]
            skip_rf = clean_values["skip_rf"]
            skip_c = clean_values["skip_c"]
            skip_j = clean_values["skip_j"]
            out_folder = clean_values["output"]
            save = clean_values["save"]
            project_name = clean_values["project_name"]

            names = [names[i] for i in range(len(files)) if files[i] != ""]
            files = [files[i] for i in range(len(files)) if files[i] != ""]

            names = [
                names[i] if names[i] != "" else files[i].split("/")[-1].split(".")[0]
                for i in range(len(files))
            ]

            # TODO: if len(files) == 1: skip_c and skip_t set to True
            # TODO:      or do something where it produces a bar graph instead??
            # TODO: also add a catch if no files are passed
            # TODO: there's a problem where RF feature importance for 2 classes produces slightly different plots
            # due to random states. Can I assign a random seed somehow?

            return (
                files,
                names,
                show,
                skip_t,
                skip_rf,
                skip_c,
                skip_j,
                out_folder,
                project_name,
                save,
            )


def calculateFeatures(files, names):
    # set up dictionary with each feature as a key and each value as an empty list for use later.
    features = {
        "classification": [],
        "name": [],
        "sequence": [],
        "length": [],
        "isoelectric point": [],
        "isoelectric point q1": [],
        "isoelectric point q2": [],
        "isoelectric point q3": [],
        "isoelectric point q4": [],
        "isoelectric point q5": [],
        "instability index": [],
        "acidic fraction": [],
        "acidic fraction q1": [],
        "acidic fraction q2": [],
        "acidic fraction q3": [],
        "acidic fraction q4": [],
        "acidic fraction q5": [],
        "basic fraction": [],
        "basic fraction q1": [],
        "basic fraction q2": [],
        "basic fraction q3": [],
        "basic fraction q4": [],
        "basic fraction q5": [],
        "molar extinction coefficient 1": [],
        "molar extinction coefficient 2": [],
        "gravy score": [],
        "gravy score q1": [],
        "gravy score q2": [],
        "gravy score q3": [],
        "gravy score q4": [],
        "gravy score q5": [],
        "flexibility score mean": [],
        "flexibility score mean q1": [],
        "flexibility score mean q2": [],
        "flexibility score mean q3": [],
        "flexibility score mean q4": [],
        "flexibility score mean q5": [],
        "peak flexibility smoothed": [],
        "peak flexibility relative position": [],
        "helix fraction": [],
        "turn fraction": [],
        "sheet fraction": [],
        "mean flexibility at ends": [],
        "molecular weight": [],
        "aromaticity index": [],
    }

    print("opening files")

    # open each files in current directory if they match the correct file type.
    count = 0
    for z in range(len(files)):
        filename = files[z]
        classification = names[z]

        with open(filename, "r") as file:
            file_string = file.read()

        print("\nOpening file: " + filename, end="\t\t")

        # convert file into a list of each fasta sequence and remove line breaks
        # zeroth index removed to remove a blank line that is somehow included and breaks things.
        fasta_list = [x.strip("\n") for x in file_string.split(">")][1:]
        print("count: " + str(len(fasta_list)) + "\n")

        # loop through each fasta name+sequence in the file
        with progressbar.ProgressBar(max_value=len(fasta_list)) as bar:
            for j in range(len(fasta_list)):
                sequence = fasta_list[j]
                count += 1
                # assign classification
                features["classification"].append(classification)

                # split the fasta sequence into a list of each line
                sequence_lines = sequence.split("\n")

                # add the "name" line to the dictionary
                features["name"].append(sequence_lines[0])

                # compress the lines after the name into a sequence string
                # and add name and length to dict
                seq = "".join(sequence_lines[1:])
                features["sequence"].append(seq)
                length = len(seq)
                features["length"].append(length)

                # cleanup the string so ProteinAnalysis can handle it
                seq = seq.upper()
                seq = (
                    seq.replace("X", "A")
                    .replace("J", "L")
                    .replace("*", "A")
                    .replace("Z", "E")
                    .replace("B", "D")
                    .replace("U", "C")
                    .replace("O", "K")
                )
                X = ProteinAnalysis(seq)

                # set up quintiles
                quintile_size = int(round(length / 5, 0))
                seq_q1, seq_q2, seq_q3, seq_q4, seq_q5 = (
                    seq[:quintile_size],
                    seq[quintile_size : quintile_size * 2],
                    seq[quintile_size * 2 : quintile_size * 3],
                    seq[quintile_size * 3 : quintile_size * 4],
                    seq[quintile_size * 4 :],
                )

                X1, X2, X3, X4, X5 = (
                    ProteinAnalysis(seq_q1),
                    ProteinAnalysis(seq_q2),
                    ProteinAnalysis(seq_q3),
                    ProteinAnalysis(seq_q4),
                    ProteinAnalysis(seq_q5),
                )

                # isoelectric point plus quintiles
                features["isoelectric point"].append(X.isoelectric_point())
                features["isoelectric point q1"].append(X1.isoelectric_point())
                features["isoelectric point q2"].append(X2.isoelectric_point())
                features["isoelectric point q3"].append(X3.isoelectric_point())
                features["isoelectric point q4"].append(X4.isoelectric_point())
                features["isoelectric point q5"].append(X5.isoelectric_point())

                # instability index
                features["instability index"].append(X.instability_index())

                # aromaticity index
                features["aromaticity index"].append(X.aromaticity())

                # molar extinction coefficients
                features["molar extinction coefficient 1"].append(
                    X.molar_extinction_coefficient()[0]
                )
                features["molar extinction coefficient 2"].append(
                    X.molar_extinction_coefficient()[1]
                )

                # gravy scores plus quintiles
                features["gravy score"].append(X.gravy())
                features["gravy score q1"].append(X1.gravy())
                features["gravy score q2"].append(X2.gravy())
                features["gravy score q3"].append(X3.gravy())
                features["gravy score q4"].append(X4.gravy())
                features["gravy score q5"].append(X5.gravy())

                # molecular weight
                features["molecular weight"].append(X.molecular_weight())

                # secondary structure fractions
                features["helix fraction"].append(X.secondary_structure_fraction()[0])
                features["turn fraction"].append(X.secondary_structure_fraction()[1])
                features["sheet fraction"].append(X.secondary_structure_fraction()[2])

                # acidic fraction plus quintiles
                features["acidic fraction"].append(
                    (seq.count("R") + seq.count("L") + seq.count("K")) / length
                )
                features["acidic fraction q1"].append(
                    (seq_q1.count("R") + seq_q1.count("L") + seq_q1.count("K"))
                    / len(seq_q1)
                )
                features["acidic fraction q2"].append(
                    (seq_q2.count("R") + seq_q2.count("L") + seq_q2.count("K"))
                    / len(seq_q2)
                )
                features["acidic fraction q3"].append(
                    (seq_q3.count("R") + seq_q3.count("L") + seq_q3.count("K"))
                    / len(seq_q3)
                )
                features["acidic fraction q4"].append(
                    (seq_q4.count("R") + seq_q4.count("L") + seq_q4.count("K"))
                    / len(seq_q4)
                )
                features["acidic fraction q5"].append(
                    (seq_q5.count("R") + seq_q5.count("L") + seq_q5.count("K"))
                    / len(seq_q5)
                )

                # basic fraction plus quintiles
                features["basic fraction"].append(
                    (seq.count("D") + seq.count("E")) / length
                )
                features["basic fraction q1"].append(
                    (seq_q1.count("D") + seq_q1.count("E")) / len(seq_q1)
                )
                features["basic fraction q2"].append(
                    (seq_q2.count("D") + seq_q2.count("E")) / len(seq_q2)
                )
                features["basic fraction q3"].append(
                    (seq_q3.count("D") + seq_q3.count("E")) / len(seq_q3)
                )
                features["basic fraction q4"].append(
                    (seq_q4.count("D") + seq_q4.count("E")) / len(seq_q4)
                )
                features["basic fraction q5"].append(
                    (seq_q5.count("D") + seq_q5.count("E")) / len(seq_q5)
                )

                # flexibility score means
                features["flexibility score mean"].append(np.mean(X.flexibility()))
                features["flexibility score mean q1"].append(np.mean(X1.flexibility()))
                features["flexibility score mean q2"].append(np.mean(X2.flexibility()))
                features["flexibility score mean q3"].append(np.mean(X3.flexibility()))
                features["flexibility score mean q4"].append(np.mean(X4.flexibility()))
                features["flexibility score mean q5"].append(np.mean(X5.flexibility()))

                # flexibility smoothed
                flex_list = X.flexibility()
                flex_smoothed = [
                    np.mean(flex_list[i : i + 5]) for i in range(len(flex_list) - 5)
                ]
                features["peak flexibility smoothed"].append(max(flex_smoothed))

                # peak flexibility relative position
                peak_index = flex_list.index(max(flex_list))
                features["peak flexibility relative position"].append(
                    peak_index / len(flex_list)
                )

                # flexibility at ends
                five_percent_range = int(round(len(flex_list) * 0.05, 0))
                start_flexibility_mean = np.mean(flex_list[:five_percent_range])
                end_flexibility_mean = np.mean(flex_list[-five_percent_range:])
                features["mean flexibility at ends"].append(
                    (start_flexibility_mean + end_flexibility_mean) / 2
                )

                bar.update(j)

    df = pd.DataFrame.from_dict(features)

    return df


def cat_cont_correlation_ratio(categories, values):
    """
   Simple function to determine the correlation ratio between a list
   of categorical values and a list of continuous values.
   Code provided by Julien.
   """
    f_cat, _ = pd.factorize(categories)
    cat_num = np.max(f_cat) + 1
    y_avg_array = np.zeros(cat_num)
    n_array = np.zeros(cat_num)
    for i in range(0, cat_num):
        cat_measures = values[np.argwhere(f_cat == i).flatten()]
        n_array[i] = len(cat_measures)
        y_avg_array[i] = np.average(cat_measures)
    y_total_avg = np.sum(np.multiply(y_avg_array, n_array)) / np.sum(n_array)
    numerator = np.sum(
        np.multiply(n_array, np.power(np.subtract(y_avg_array, y_total_avg), 2))
    )
    denominator = np.sum(np.power(np.subtract(values, y_total_avg), 2))
    if numerator == 0:
        eta = 0.0
    else:
        eta = np.sqrt(numerator / denominator)
    return eta


def correlation_bar_plots(df, outfile, class_list, show):
    # get a list of the name of each feature if the first value in the dataframe in the target column is a number
    features_list = [
        x
        for x in df.columns.to_list()
        if type(df[x][0]) == np.float64 or type(df[x][0]) == np.int64
    ]
    # sort it alphabetically ignoring case
    features_list = sorted(features_list, key=str.lower)

    # set up an empty array for correlation values. Rows = classes, Columns = features
    correlation_array = np.zeros((len(class_list), len(features_list)))

    for i in range(len(class_list)):
        classification = class_list[i]

        # figure out which rows are the correct classification
        is_target_class = [
            "a" if x == classification else "b" for x in df["classification"]
        ]

        # loop through each feature from the list
        for ii in range(len(features_list)):
            feature = features_list[ii]
            # run the correlation function between it and the classification column
            correlation_array[i][ii] = round(
                cat_cont_correlation_ratio(is_target_class, df[feature]), 5
            )

    correlation_plot = go.Figure(
        data=go.Heatmap(
            z=correlation_array,
            x=features_list,
            y=class_list,
            hoverongaps=False,
            colorscale="Plasma",
        )
    )

    correlation_plot.update_layout(
        xaxis_title="Features",
        yaxis_title="Phage Protein Class",
        title_text=f"Correlation Ratios",
        font=dict(size=12),
    )
    correlation_plot.update_xaxes(tickangle=45, tickfont=dict(size=12))

    correlation_plot.write_html(
        file=f"{outfile}/correlation heatplot.html", include_plotlyjs=True
    )

    if show:
        correlation_plot.show()

    return


def t_value_bar_plots(df, outfile, class_list, show):
    # get a list of the name of each feature if the first value in the dataframe in the target column is a number
    features_list = [
        x
        for x in df.columns.to_list()
        if type(df[x][0]) == np.float64 or type(df[x][0]) == np.int64
    ]
    # sort it alphabetically ignoring case
    features_list = sorted(features_list, key=str.lower)

    # set up an empty array for correlation values. Rows = classes, Columns = features
    tval_array = np.zeros((len(class_list), len(features_list)))

    for i in range(len(class_list)):
        classification = class_list[i]

        # figure out which rows are the correct classification
        is_target_class = [
            1 if x == classification else 0 for x in df["classification"]
        ]

        # loop through each column in the data frame and check if the first row value is a number of some kind
        for ii in range(len(features_list)):
            feature = features_list[ii]
            if type(df[feature][0]) == np.float64 or type(df[feature][0]) == np.int64:
                # if it is, get the t-value. Following code was provided by Julien.
                predictor = statsmodels.api.add_constant(df[feature].to_numpy())

                logistic_regression_model = statsmodels.api.Logit(
                    is_target_class, predictor
                )
                logistic_regression_fitted = logistic_regression_model.fit(disp=False)

                t_value = round(logistic_regression_fitted.tvalues[1], 4)
                tval_array[i][ii] = abs(t_value)

    t_val_plot = go.Figure(
        data=go.Heatmap(
            z=tval_array,
            x=features_list,
            y=class_list,
            hoverongaps=False,
            colorscale="Plasma",
        )
    )

    t_val_plot.update_layout(
        xaxis_title="Features",
        yaxis_title="Phage Protein Class",
        title_text=f"T-Values",
        font=dict(size=12),
    )
    t_val_plot.update_xaxes(tickangle=45, tickfont=dict(size=12))

    t_val_plot.write_html(
        file=f"{outfile}/t-value heatplot.html", include_plotlyjs=True
    )

    if show:
        t_val_plot.show()

    return


def feature_importance_bar_plot(df, outfile, class_list, show):
    for classification in class_list:

        # figure out which rows are the correct classification, convert into data frame
        is_target_class = [
            True if x == classification else False for x in df["classification"]
        ]
        is_target_class = pd.DataFrame(is_target_class, columns=["classification"])

        # setup second data frame that drops all non-testable columns
        df_rf = df
        for feature in df_rf.columns.to_list():
            if (
                type(df_rf[feature][0]) != np.float64
                and type(df_rf[feature][0]) != np.int64
            ):
                df_rf = df_rf.drop(feature, axis=1)

        # convert to numpy and flatten the is_target_class array to make RandomForestRegressor() happy
        is_target_class_array = np.ravel(is_target_class.to_numpy(), order="C")
        df_rf_array = df_rf.to_numpy()

        # get feature importance for each column in the data frame
        rf = RandomForestRegressor()
        rf.fit(df_rf_array, is_target_class_array)
        feature_importance = rf.feature_importances_

        # set x and y values for plotting
        x_axis = list(df_rf.columns)
        y_axis = list(feature_importance)

        # sort numerically using zip sort and then return data to list format
        y_axis, x_axis = zip(*sorted(zip(y_axis, x_axis)))
        x_axis, y_axis = list(x_axis), list(y_axis)

        # set up correlation plot with layout options
        rf_plot = go.Figure(
            [
                go.Bar(
                    x=x_axis, y=y_axis, marker={"color": y_axis, "colorscale": "Blugrn"}
                )
            ]
        )
        rf_plot.update_layout(
            xaxis_title="Features",
            yaxis_title="RF Feature Importance Metric",
            title_text=f"Random Forest Feature Importance - {classification}",
            font=dict(size=12),
        )
        rf_plot.update_xaxes(tickangle=45, tickfont=dict(size=10))

        rf_plot.write_html(
            file=f"{outfile}/RF feature importance {classification}.html",
            include_plotlyjs=True,
        )
        if show:
            rf_plot.show()
    return


def joy_plot(df, outfile, class_list, show):
    for feature in df.columns.to_list():
        if type(df[feature][0]) == np.float64 or type(df[feature][0]) == np.int64:
            df_to_plot = df[["classification", feature]]
            df_to_plot = df_to_plot[df["classification"].isin(class_list)]
            joy_violin = px.violin(
                df_to_plot,
                x=feature,
                color="classification",
                violinmode="overlay",
                points=False,
                orientation="h",
                color_discrete_sequence=px.colors.qualitative.Pastel1,
            )
            joy_violin.update_layout(height=400)
            joy_violin.update_traces(width=0.9, points=False)
            joy_violin.update_yaxes(range=[0, 1])

            joy_violin.update_layout(title=f"Joy Plot of {feature}")
            joy_violin.update_layout(xaxis_showgrid=False, xaxis_zeroline=False)

            joy_violin.write_html(
                file=f"{outfile}/joy plot {feature}.html", include_plotlyjs="cdn"
            )
            if show:
                joy_violin.show()
    return


def graphDriver(df, outfile, classes, show, skip_t, skip_rf, skip_c, skip_j):
    # open dataframe from csv
    df = df.dropna(axis=0).reset_index()

    if "index" in df.columns.to_list():
        df.drop("index", axis=1, inplace=True)

    # generate joy plots for each feature
    if not skip_j:
        print("Generating joy plots")
        joy_plot(df, outfile, classes, show)
    else:
        print("Skipping joy plots")

    # generate stand-alone correlation bar plots
    if not skip_c:
        print("Generating correlation plots")
        correlation_bar_plots(df=df, outfile=outfile, class_list=classes, show=show)
    else:
        print("Skipping correlation plots")

    # generate stand-alone t-value bar plots
    if not skip_t:
        print("Generating t-value plots")
        t_value_bar_plots(df=df, outfile=outfile, class_list=classes, show=show)
    else:
        print("Skipping t-value plots")

    # generate RF feature importance plots
    if not skip_rf:
        print("Generating feature importance plots. This may take some time...")
        feature_importance_bar_plot(
            df=df, outfile=outfile, class_list=classes, show=show
        )
    else:
        print("Skipping feature importance plots")

    print("Success!")


def main():
    (
        files,
        names,
        show,
        skip_t,
        skip_rf,
        skip_c,
        skip_j,
        out_dir,
        project_name,
        save,
    ) = simpleGUI()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    project_dir = out_dir + "/" + project_name + "/"
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    logfile = project_dir + f"log_{project_name}.txt"
    with open(logfile, "w") as log:
        now = time.strftime("%Y,%m,%d,%H,%M,%S").split(",")
        log.write(
            f"Project Name: {project_name}\n"
            f"Project Date: {now[1]} / {now[2]} / {now[0]}\n\n"
        )

    sys.stdout = open(logfile, "a")

    bad_files = fileQC(files)

    if len(bad_files) > 0:
        bad_files_str = "\n".join(bad_files)
        sys.exit(
            print(
                f"Errors were detected in the following files and the run was aborted:\n{bad_files_str}"
            )
        )

    if len(files) == 1:
        print(
            "Only 1 file provided. Skipping t-value plot, rf plot, and correlation plot.\n"
        )
        skip_t, skip_rf, skip_c = True, True, True

    df = calculateFeatures(files, names)

    if save:
        print("Saving File")
        if not project_name.endswith(".csv"):
            df.to_csv(project_dir + project_name + ".csv", index=False, header=True)
        else:
            df.to_csv(project_dir + project_name, index=False, header=True)

    graphDriver(df, project_dir, names, show, skip_t, skip_rf, skip_c, skip_j)

    sys.stdout.flush()


if __name__ == "__main__":
    main()
