import pandas as pd
import plotnine as p9
import numpy as np
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef

pd.options.mode.chained_assignment = None  # default='warn'


def return_output_analysis(true, pred):

    # This method writes scikit-learn's classification report and the MCC to "classification_scores.txt"

    print(classification_report(true, pred, digits=4))
    print(f"Matthews Correlation Coefficient: {np.round(matthews_corrcoef(true, pred), 4)}")

    output_file = "classification_scores.txt"
    output = open(output_file, "w")
    output.write(classification_report(true, pred, digits=4))
    output.write(f"Matthews Correlation Coefficient: {np.round(matthews_corrcoef(true, pred), 4)}")
    output.close()
    print(f"Classification scores are saved in {output_file}")


def plot_confusion_matrix(true, pred):

    # This method plots a normalized confusion matrix to a file called "confusion_matrix.png"

    conf_mat = confusion_matrix(true, pred)

    conf = []
    for row in conf_mat:
        conf.append(row / sum(row))

    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]

    conf = pd.DataFrame(conf, index=rna_types, columns=rna_types)
    conf["rna_type"] = conf.index  # add index column for melting
    rna_types = ["lncRNA", "miRNA", "rRNA", "snRNA", "snoRNA", "tRNA"]  # desired order of labels
    conf_melt = pd.melt(conf, id_vars="rna_type", ignore_index=True, value_name="Prediction",
                        var_name="Predicted Class")  # melt into df with 3 columns
    conf_melt["Predicted Class"] = pd.Categorical(conf_melt["Predicted Class"],
                                                  categories=rna_types,
                                                  ordered=True)

    text_color = np.array(['black'] * len(conf_melt))
    # change color according to box color
    text_color[conf_melt['Prediction'] > 0.5] = 'white'
    text_color[conf_melt['Prediction'] < 0.005] = "#FFF7EC"

    # This creates the labels displayed in the conf_matrix
    # For values below the rounding threshold (<0.005) the box is left empty to reduce clutter
    conf_melt["pred_string"] = ""
    for i in range(0, len(conf_melt)):
        if conf_melt.Prediction[i] >= 0.005:
            conf_melt.pred_string[i] = str(np.round(conf_melt.Prediction[i], 2))

    plot = (p9.ggplot(conf_melt, p9.aes('Predicted Class', 'rna_type', fill="Prediction"))
            + p9.geom_tile(p9.aes(width=.95, height=.95))
            + p9.geom_text(p9.aes(label='pred_string'), size=8, color=text_color, fontweight="bold") # Add labels
            + p9.scale_y_discrete(limits=list(reversed(rna_types))) # Change order according to preference
            + p9.scale_fill_distiller(type="seq", palette="OrRd", direction=1, guide=False) # Change color to gradient
            + p9.xlab("Predicted Label")
            + p9.ylab("True Label")
            + p9.labs(title="Confusion Matrix")
            + p9.theme(
                figure_size=(3, 3),
                panel_background=p9.element_rect(fill='white'),
                axis_text_x=p9.element_text(colour="black", size=8, rotation=45, ha="right"),
                axis_text_y=p9.element_text(colour="black", size=8, va="center"),
                axis_title_x=p9.element_text(colour="black", size=10),
                axis_title_y=p9.element_text(colour="black", size=10),
                plot_title=p9.element_blank()
                )
            )
    plot.save(filename="confusion_matrix.png", dpi=300)
    print("The confusion matrix was saved in results/confusion_matrix.png")
