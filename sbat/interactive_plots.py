import os
import shutil

import pandas as pd
from bokeh.models import HoverTool, DatePicker, FileInput, Button, ColumnDataSource, RangeSlider, DataTable, \
    TableColumn, NumericInput, StringFormatter, NumberFormatter
from bokeh.plotting import figure
from bokeh.server.server import Server
from tornado.ioloop import IOLoop
from bokeh.layouts import column
from bokeh.models import Button, FileInput, NumericInput, RangeSlider

PLOTTER = None
import utils


class IPlotter:
    def __init__(self, name, size, data, margin=5):
        self.size = size
        self.margin = margin
        self.filename = name
        self.stats = pd.read_csv("np_sb_analysis_example.csv")
        self.df = data
        self.df["more_freq_count"] = self.df.apply(lambda row: utils.select_more_frequent(row),
                                                                         axis=1)
        self.df["selected"] = False
        self.df["index"] = [i for i in range(len(self.df))]
        self.max_freq = self.df["more_freq_count"].max()
        self.ds = ColumnDataSource().data = self.df
        self.plot_ds = ColumnDataSource(self.df)

        self.kmer_vs_sb_df = self.df.copy()
        self.kmer_vs_sb_ds = ColumnDataSource(self.kmer_vs_sb_df)
        self.nanopore = False
        self.data_table = self.create_data_table()
        self.gc_data = None

    def create_gc_plot(self):
        plot = figure(
            plot_height=500,
            plot_width=1000,
            title="Strand Bias Levels in Relation to GC Percentage",
            tools="crosshair, pan,reset, save,wheel_zoom, box_select, "
                  "poly_select, tap, box_zoom",
        )
        print("pased plot creation")

        tooltips = [
            ("Average GC Content (%)", "@x"),
            ("Average Strand Bias", "@y"),
        ]

        plot.add_tools(
            HoverTool(
                tooltips=tooltips,
            )
        )

        plot.inverted_triangle(self.gc_data.upper_gc,
                      self.gc_data.upper_biases,
                      color="red",
                      fill_color="red",
                      size=10,
                      legend_label="GC Content vs Strand Bias in Top {}% of SB Score".format(self.margin))
        print("passed triangle 1")
        plot.triangle(self.gc_data.lower_gc,
                      self.gc_data.lower_biases,
                      color="green",
                      fill_color="green",
                      size=10,
                      legend_label="GC Content vs Strand Bias in Bottom {}% of SB Score".format(self.margin))
        print("passed triangle 2")
        plot.legend.location = "top_left"
        #self.gc_plot = plot
        return plot

    def create_lineplot(self, bins=True, k=None):
        if bins:
            x_axis = "Bin"
            title = 'Statistics of Strand Bias among Bins for K={}'.format(k)
        else:
            x_axis = "K"
            title = 'Mean and Median of Strand Bias for Different Sizes of K'
        plot = figure(
            plot_height=500,
            plot_width=800,
            title=title,
            tools="crosshair, pan,reset, save,wheel_zoom, box_select, "
                  "poly_select, tap, box_zoom",
        )

        ds = ColumnDataSource(self.df)
        mean = plot.line(self.stats[x_axis.lower()], self.stats["bias_mean"], line_color="green", legend_label="Mean", line_width=3)
        median = plot.line(self.stats[x_axis.lower()], self.stats["bias_median"], line_color="blue", legend_label="Median", line_width=3)
        modus = plot.line(self.stats[x_axis.lower()], self.stats["bias_modus"], line_color="pink", legend_label="Mode", line_width=3)
        perc_5 = plot.line(self.stats[x_axis.lower()], self.stats["percentile_5"], line_color="orange", legend_label="5th Percentile", line_width=3)
        perc_95 = plot.line(self.stats[x_axis.lower()], self.stats["percentile_95"], line_color="red", legend_label="95th Percentile", line_width=3)

        plot.add_tools(HoverTool(tooltips="Median Strand Bias: @y {}: @x".format(x_axis), renderers=[median], mode="vline"))
        plot.add_tools(HoverTool(tooltips="Mean Strand Bias: @y {}: @x".format(x_axis), renderers=[mean], mode="vline"))
        plot.add_tools(HoverTool(tooltips="Mode Strand Bias: @y {}: @x".format(x_axis), renderers=[modus], mode="vline"))
        plot.add_tools(HoverTool(tooltips="5th Percentile of Strand Bias: @y {}: @x".format(x_axis), renderers=[perc_5], mode="vline"))
        plot.add_tools(HoverTool(tooltips="95th Percentile of Strand Bias: @y {}: @x".format(x_axis), renderers=[perc_95], mode="vline"))

        plot.legend.click_policy = "hide"
        plot.xaxis.axis_label = "Bins"
        plot.yaxis.axis_label = "Strand Bias [%]"
        plot.xaxis.axis_label_text_font_size = "15pt"
        plot.yaxis.axis_label_text_font_size = "15pt"
        plot.title.text_font_size = "15pt"
        return plot


    def create_kmers_vs_bias_plot(self):
        """
        Creates plot figure and adds tooltips to it

        :return:
        """
        plot = figure(
            plot_height=800,
            plot_width=1400,
            name="kmers_vs_sb",
            title="Kmers vs Strand Bias",
            tools="crosshair, pan,reset, save,wheel_zoom, box_select, "
            "poly_select, tap, box_zoom, lasso_select",
        )

        tooltips = [
            ("K-mer", "@seq"),
            ("Frequency", "@seq_count"),
            ("Complement", "@rev_complement"),
            ("Frequency", "@rev_complement_count"),
            ("Strand Bias (%)", "@{strand_bias_%}"),
            ("GC Content (%)", "@{GC_%}"),
        ]


        plot.add_tools(
            HoverTool(
                tooltips=tooltips,
            )
        )

        # get more frequent out of k-mer and its rev. complement
        self.kmer_vs_sb_df = self.kmer_vs_sb_df.sort_values(by=['more_freq_count'], ascending=False, ignore_index=True)
        self.kmer_vs_sb_df["order_index"] = range(1, len(self.kmer_vs_sb_df) + 1)
        self.kmer_vs_sb_ds = ColumnDataSource(self.kmer_vs_sb_df)

        plot.circle(
            "order_index",
            "strand_bias_%",
            source=self.kmer_vs_sb_ds,
            name="Strand bias of k-mers of size {} in relation with frequency".format(self.size),
            legend_label="Strand bias of k-mers of size {} in relation with frequency".format(self.size),
            color="green",
        )

        plot.yaxis.axis_label = "Strand Bias [%]"
        plot.xaxis.axis_label_text_font_size = "15pt"
        plot.yaxis.axis_label_text_font_size = "15pt"
        plot.title.text_font_size = "15pt"
        self.kmers_vs_sb = plot
        return plot


    def create_data_table(self):
        columns = [
            TableColumn(
                field="seq",
                title="K-mer",
                formatter=StringFormatter(),
                width=200
            ),
            TableColumn(
                field="seq_count",
                title="K-mer Frequency",
                formatter=NumberFormatter(format="0,0"),
                width=400
            ),
            TableColumn(
                field="rev_complement",
                title="Reverse Complement",
                formatter=StringFormatter(),
            ),
            TableColumn(
                field="rev_complement_count",
                title="Complement Frequency",
                formatter=NumberFormatter(format="0,0"),
            ),
            TableColumn(
                field="strand_bias_%",
                title="Strand Bias (%)",
                formatter=NumberFormatter(format="0[.]0000"),
            ),
            TableColumn(
                field="GC_%",
                title="GC Content (%)",
                formatter=NumberFormatter(format="0[.]00"),
            )
        ]

        return DataTable(columns=columns, width=900)

    def download_selected(self):
        """
        Is triggered by pushing download button. Saves dataframe to ./labeled folder
        under loaded file filename + "-labeled.csv" extension. Filters dataframe
        according to "days to download" field.

        :return:
        """

        download_df = self.kmer_vs_sb_df[self.kmer_vs_sb_df["selected"]]
        if not os.path.isdir("selected"):
            os.makedirs("selected")

        filename = self.filename + "-labeled.csv"
        download_df.to_csv(
            os.path.join("selected/", filename),
            index=False,
        )
        print("Downloading dataset.")

    def update_selected_in_df(self, attrname, old, new):
        """
        Called upon selecting datapoints from marked anomaly on graph. Function updates
        selected column in dataframe.

        :param attrname:
        :param old:
        :param new:
        :return:
        """
        self.kmer_vs_sb_df["selected"] = False
        for index in new:
            df_index = self.kmer_vs_sb_ds.data["order_index"][index]
            self.kmer_vs_sb_df.at[df_index, "selected"] = True
        selected_df = self.kmer_vs_sb_df[self.kmer_vs_sb_df["selected"]]
        self.data_table.source.data = selected_df

    def clear_selected(self):
        """
        Is called upon pressing clear selected button. Clears table and graph with
        selected data points. This function does not modify working dataframe.
        :return:
        """
        self.kmer_vs_sb_ds.selected.indices = []

        self.data_table.source.selected.indices = []
        self.data_table.source.data = pd.DataFrame()
        print("Clearing selected in graph and datatable.")


def remove_glyphs(plot, plot_names=None):
    """
    Remove glyphs from plot by given names

    :param plot:
    :param plot_names:
    :return:
    """
    for name in plot_names:
        try:
            plot.renderers.remove(name)
        except ValueError:
            pass





days_for_visual_inspection = NumericInput(
    value=0,
    title="Days for visual inspection:",
    low=0,
    high=10,
)
days_to_label_input = NumericInput(
    value=5,
    title="Days to label:",
    low=0,
    high=10,
)
file_button = FileInput(accept=".csv")
download_button = Button(label="Download selected", button_type="success")
clear_selected_button = Button(label="Clear selected", button_type="primary")


def analyze_bias(analysis, jf, nano, input_files):
    """
    Function executing the core of the application.
    """
    for file in input_files:  # run SBAT for each input file
        analysis.set_file(file)
        analysis.init_analysis()
        print("input: " + file)
        dfs = []
        detect_nano = nano is not None and utils.check_if_nanopore(file)
        run_nano = detect_nano
        for k in range(analysis.start_k, analysis.end_k + 1):
            if jf is not None:
                print("running computation and analysis for K=" + str(k))
                jf_output = jf.run_jellyfish(file, k)

                if run_nano:
                    print('running nanopore analysis...')
                    nano.nanopore_analysis(file)
                    run_nano = False  # run analysis only once for each input file
            else:
                print("jellyfish disabled, running only analysis...")
                jf_output = file
            df = analysis.jellyfish_to_dataframe(jf_output, k)  # convert jellyfish results to DataFrame
            dfs.append(df)
            if df is not None:
                analysis.plot_kmers_vs_bias(df, k)
                if analysis.interactive and analysis.plotter is None:
                    analysis.plotter = IPlotter(analysis.filename, k, df)
                    global PLOTTER
                    PLOTTER = analysis.plotter
                if not detect_nano:
                    # if not nanopore, create simplified version of this stats,
                    # otherwise it is created as part of nanopore analysis
                    analysis.track_most_common_kmer_change_freq([df], k)
        if analysis.interactive:
            analysis.plotter.gc_data = analysis.calculate_gc_plot_data(dfs)
            analysis.plot_basic_stats_lineplot(analysis.filename, analysis.sb_analysis_file)
            analysis.plot_conf_interval_graph(dfs, start_index=analysis.start_k)
            analysis.plot_gc_from_dataframe(analysis.plotter.gc_data)
    if not analysis.keep_computations:
        shutil.rmtree(analysis.dump_dir)
        if jf is not None:
            shutil.rmtree(jf.jf_dir)




def interactive_plot(doc):
    if PLOTTER is None:
        return

    clear_selected_button.on_click(PLOTTER.clear_selected)
    range_slider_sb = RangeSlider(title="Filter K-mers by Strand Bias", start=0, end=100, value=(0, 100), step=1)
    range_slider_frequency = RangeSlider(title="Filter K-mers by Frequency", start=0, end=PLOTTER.max_freq, value=(0, PLOTTER.max_freq), step=max(1000, PLOTTER.max_freq//100))
    # days_to_label_input.on_change("value", plotter.update_data_to_label)

    def range_refresh(attr, old, new):
        PLOTTER.kmer_vs_sb_ds.data = PLOTTER.kmer_vs_sb_df[
            (PLOTTER.kmer_vs_sb_df["strand_bias_%"] >= range_slider_sb.value[0]) &
            (PLOTTER.kmer_vs_sb_df["strand_bias_%"] <= range_slider_sb.value[1]) &
            (PLOTTER.kmer_vs_sb_df["more_freq_count"] >= range_slider_frequency.value[0]) &
            (PLOTTER.kmer_vs_sb_df["more_freq_count"] <= range_slider_frequency.value[1])
            ].sort_values(by="index")

    range_slider_sb.on_change("value", range_refresh)#.on_change("value", plotter.refresh)
    range_slider_frequency.on_change("value", range_refresh)
    download_button.on_click(PLOTTER.download_selected)

    lineplot = PLOTTER.create_lineplot()
    kmers_vs_sb = PLOTTER.create_kmers_vs_bias_plot()
    print(PLOTTER.gc_data, PLOTTER.margin)
    gc_plot = PLOTTER.create_gc_plot()

    PLOTTER.kmer_vs_sb_ds.selected.on_change(
          "indices", PLOTTER.update_selected_in_df
    )

    buttons = column(clear_selected_button,range_slider_sb, range_slider_frequency)
    download = column(download_button)
    doc.add_root(column(children=[lineplot, gc_plot, buttons, kmers_vs_sb, download, PLOTTER.data_table], width=1400))
    return doc

def run_server():
    server = Server({'/bkapp': interactive_plot}, io_loop=IOLoop())

    server.start()
    from bokeh.util.browser import view
    #
    server.io_loop.add_callback(view, "http://localhost:5006/")
    server.io_loop.start()
    server.io_loop.close()