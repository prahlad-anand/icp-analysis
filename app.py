"""
Interactive Dash dashboard for immune cell population analysis.

Run:
    1. python app.py
    2. Open http://127.0.0.1:8050 in a browser.
"""

import sqlite3
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, callback, dash_table, dcc, html
from scipy import stats


DB_PATH = "cell_data.db"
POPULATIONS = ["b_cell", "cd8_t_cell", 'cd4_t_cell', "nk_cell", "monocyte"]
LABELS2TXT = {
    "b_cell":      "B Cell",
    "cd8_t_cell":  "CD8 T Cell",
    "cd4_t_cell":  "CD4 T Cell",
    "nk_cell":     "NK Cell",
    "monocyte":    "Monocyte",
}
RESPONSE_COLORS = {"yes": "#00cc00", "no": "#cc0000"}
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-Responder"}


def query(sql):
    conn = sqlite3.connect(DB_PATH)
    df = pd.read_sql_query(sql, conn)
    conn.close()

    return df


"""
PART 2: Frequency data
"""
def build_frequency_df():
    """
    Return a long-format DataFrame with one row per (sample, population) pair.
    """
    raw = query(
        """
        SELECT
            s.sample_id   AS sample,
            su.project_id AS project,
            su.condition,
            su.treatment,
            s.sample_type,
            s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        """
    )
    raw["total_count"] = raw[POPULATIONS].sum(axis=1)

    long = raw.melt(
        id_vars=["sample", "project", "condition", "treatment", "sample_type", "total_count"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["percentage"] = (long["count"] / long["total_count"] * 100).round(2)

    return long


"""
PART 3: Statistical Analysis
"""
def build_stat_data():
    """
    Returns (plot_df, stat_df) for the melanoma & miraclib & PBMC cohort.
    plot_df – long format for boxplot (sample_id, response, population, percentage)
    stat_df – one row per population with t-test results
    """
    raw = query(
        """
        SELECT
            s.sample_id,
            su.response,
            s.b_cell, s.cd8_t_cell, s.cd4_t_cell, s.nk_cell, s.monocyte
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        WHERE su.condition = 'melanoma'
        AND su.treatment = 'miraclib'
        AND s.sample_type = 'PBMC'
        AND su.response IS NOT NULL
        """
    )

    raw["total"] = raw[POPULATIONS].sum(axis=1)
    pct_cols = []
    for p in POPULATIONS:
        col = f"{p}_pct"
        raw[col] = (raw[p] / raw["total"] * 100).round(2)
        pct_cols.append(col)

    plot_df = raw.melt(
        id_vars=["sample_id", "response"],
        value_vars=pct_cols,
        var_name="population",
        value_name="percentage",
    )
    plot_df["population"] = plot_df["population"].str.replace("_pct", "", regex=False)
    plot_df["pop_label"]  = plot_df["population"].map(LABELS2TXT)
    plot_df["resp_label"] = plot_df["response"].map(RESPONSE_LABELS)

    # Welch's t-test per population
    rows = []
    for p in POPULATIONS:
        col = f"{p}_pct"
        yes_vals = raw.loc[raw["response"] == "yes", col].values
        no_vals  = raw.loc[raw["response"] == "no",  col].values

        t_stat, pval = stats.ttest_ind(yes_vals, no_vals, equal_var=False)

        rows.append(
            {
                "Population": LABELS2TXT[p],
                "Responders (n)": int(len(yes_vals)),
                "Non-Responders (n)": int(len(no_vals)),
                "Mean – Responders (%)": round(float(np.mean(yes_vals)), 2),
                "Mean – Non-Responders (%)": round(float(np.mean(no_vals)),  2),
                "t-statistic": round(float(t_stat), 3),
                "p-value": round(float(pval), 4),
                "Significant (p < 0.05)":   "Yes" if pval < 0.05 else "No",
            }
        )

    return plot_df, pd.DataFrame(rows)


"""
PART 4: Subset Analysis
"""
def build_subset_data():
    """
    Query: melanoma & miraclib & PBMC & time_from_treatment_start = 0

    Returns a tuple of datadrames:
    (full_df, samples_per_project, subjects_by_response, subjects_by_sex,n_samples, n_subjects)
    """
    df = query(
        """
        SELECT
            s.sample_id,
            su.subject_id,
            su.project_id,
            su.response,
            su.sex
        FROM samples s
        JOIN subjects su ON s.subject_id = su.subject_id
        WHERE su.condition = 'melanoma'
        AND su.treatment = 'miraclib'
        AND s.sample_type = 'PBMC'
        AND s.time_from_treatment_start = 0
        """
    )

    samples_per_project = (
        df.groupby("project_id", sort=True)
        .size()
        .reset_index(name="Sample Count")
        .rename(columns={"project_id": "Project"})
    )

    subj = df.drop_duplicates("subject_id")

    subjects_by_response = (
        subj.groupby("response", sort=True)
        .size()
        .reset_index(name="Subject Count")
        .rename(columns={"response": "Response"})
    )
    subjects_by_response["Response"] = subjects_by_response["Response"].map(
        {"yes": "Responder", "no": "Non-Responder"}
    )

    subjects_by_sex = (
        subj.groupby("sex", sort=True)
        .size()
        .reset_index(name="Subject Count")
        .rename(columns={"sex": "Sex"})
    )
    subjects_by_sex["Sex"] = subjects_by_sex["Sex"].map({"M": "Male", "F": "Female"})

    return df, samples_per_project, subjects_by_response, subjects_by_sex, len(df), len(subj)

print("Loading data…")
freq_df = build_frequency_df()
part3_plot_df, part3_stat_df = build_stat_data()
_, proj_df, resp_df, sex_df, n_samples_p4, n_subjects_p4 = build_subset_data()
print("Data loaded.")

# Filters for Part 2
_projects      = sorted(freq_df["project"].unique())
_conditions    = sorted(freq_df["condition"].unique())
_treatments    = sorted(freq_df["treatment"].unique())
_sample_types  = sorted(freq_df["sample_type"].unique())


# ── Shared style helpers ───────────────────────────────────────────────────────

CARD = {
    "backgroundColor": "#ffffff",
    "padding": "24px",
    "boxShadow": "0 2px 8px rgba(0,0,0,0.08)",
    "marginBottom": "24px",
}

TABLE_HEADER = {"backgroundColor": "#ff0000", "color": "white", "fontWeight": "bold"}
TABLE_CELL   = {
    "fontFamily": "'Segoe UI', sans-serif",
    "fontSize": "13px",
    "padding": "8px 12px",
    "border": "1px solid #e0e0e0",
}
TABLE_ODD    = [{"if": {"row_index": "odd"}, "backgroundColor": "#ffffff"}]


def dropdown_opts(values):
    return [{"label": "All", "value": "__all__"}] + [{"label": v, "value": v} for v in values]


def mini_table(df: pd.DataFrame, table_id: str):
    return dash_table.DataTable(
        id=table_id,
        data=df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],
        style_header=TABLE_HEADER,
        style_cell={**TABLE_CELL, "textAlign": "center"},
        style_data_conditional=TABLE_ODD,
    )

def make_boxplot(df: pd.DataFrame) -> go.Figure:
    fig = px.box(
        df,
        x="pop_label",
        y="percentage",
        color="resp_label",
        color_discrete_map={
            RESPONSE_LABELS["yes"]: RESPONSE_COLORS["yes"],
            RESPONSE_LABELS["no"]:  RESPONSE_COLORS["no"],
        },
        points="all",
        labels={
            "pop_label":   "Cell Population",
            "percentage":  "Relative Frequency (%)",
            "resp_label":  "Response",
        },
        title="Cell Population Relative Frequencies: Responders vs Non-Responders",
        category_orders={
            "pop_label":  [LABELS2TXT[p] for p in POPULATIONS],
            "resp_label": [RESPONSE_LABELS["yes"], RESPONSE_LABELS["no"]],
        },
    )
    fig.update_traces(marker_size=4, jitter=0.3)
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": "'Segoe UI', sans-serif", "size": 13},
        legend_title_text="Response",
        xaxis={"showgrid": False},
        yaxis={"gridcolor": "#ffffff", "zeroline": False},
    )

    return fig


"""
App
"""
app = Dash(__name__)
app.title = "Loblaw Bio: Immune Cell Analysis"

app.layout = html.Div(
    [
        # Header
        html.Div(
            [
                html.H1(
                    "Loblaw Bio",
                    style={"color": "#ff0000", "marginBottom": "4px", "fontSize": "28px"},
                ),
                html.P(
                    "Immune Cell Population Analysis",
                    style={"color": "#cc0000", "margin": 0, "fontSize": "14px"},
                ),
            ],
            style={
                "background-color": "#000000",
                "padding": "24px 40px",
            },
        ),

        # Tabs for each part
        html.Div(
            dcc.Tabs(
                [
                    # Part 2
                    dcc.Tab(
                        label="Part 2: Frequency Overview",
                        children=[
                            html.Div(
                                [
                                    html.H3("Initial Analysis - Data Overview", style={"color": "#cc0000"},),
                                    html.P(
                                        "Cell Population Relative Frequencies"
                                    ),

                                    # Filters
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Label("Project", style={"fontWeight": "600"}),
                                                    dcc.Dropdown(
                                                        id="f-project",
                                                        options=dropdown_opts(_projects),
                                                        value="__all__",
                                                        clearable=False,
                                                    ),
                                                ],
                                                style={"flex": 1, "minWidth": "130px", "marginRight": "12px"},
                                            ),
                                            html.Div(
                                                [
                                                    html.Label("Condition", style={"fontWeight": "600"}),
                                                    dcc.Dropdown(
                                                        id="f-condition",
                                                        options=dropdown_opts(_conditions),
                                                        value="__all__",
                                                        clearable=False,
                                                    ),
                                                ],
                                                style={"flex": 1, "minWidth": "130px", "marginRight": "12px"},
                                            ),
                                            html.Div(
                                                [
                                                    html.Label("Treatment", style={"fontWeight": "600"}),
                                                    dcc.Dropdown(
                                                        id="f-treatment",
                                                        options=dropdown_opts(_treatments),
                                                        value="__all__",
                                                        clearable=False,
                                                    ),
                                                ],
                                                style={"flex": 1, "minWidth": "130px", "marginRight": "12px"},
                                            ),
                                            html.Div(
                                                [
                                                    html.Label("Sample Type", style={"fontWeight": "600"}),
                                                    dcc.Dropdown(
                                                        id="f-sample-type",
                                                        options=dropdown_opts(_sample_types),
                                                        value="__all__",
                                                        clearable=False,
                                                    ),
                                                ],
                                                style={"flex": 1, "minWidth": "130px", "marginRight": "12px"},
                                            ),
                                            html.Div(
                                                [
                                                    html.Label("Population", style={"fontWeight": "600"}),
                                                    dcc.Dropdown(
                                                        id="f-population",
                                                        options=dropdown_opts(POPULATIONS),
                                                        value="__all__",
                                                        clearable=False,
                                                    ),
                                                ],
                                                style={"flex": 1, "minWidth": "130px"},
                                            ),
                                        ],
                                        style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "marginBottom": "16px"},
                                    ),

                                    html.Div(id="freq-info", style={"color": "#000000", "marginBottom": "8px", "fontSize": "13px"}),

                                    dash_table.DataTable(
                                        id="freq-table",
                                        columns=[
                                            {"name": "Sample",          "id": "sample"},
                                            {"name": "Total Count",     "id": "total_count",  "type": "numeric"},
                                            {"name": "Population",      "id": "population"},
                                            {"name": "Count",           "id": "count",        "type": "numeric"},
                                            {"name": "Percentage (%)",  "id": "percentage",   "type": "numeric",
                                             "format": {"specifier": ".2f"}},
                                        ],
                                        page_size=20,
                                        sort_action="native",
                                        style_table={"overflowX": "auto"},
                                        style_header=TABLE_HEADER,
                                        style_cell={**TABLE_CELL, "fontFamily": "monospace"},
                                        style_data_conditional=TABLE_ODD,
                                    ),
                                ],
                                style=CARD,
                            )
                        ],
                    ),

                    # Part 3
                    dcc.Tab(
                        label="Part 3: Responders vs Non-Responders",
                        children=[
                            html.Div(
                                [
                                    html.H3("Statistical Comparison: Responders vs Non-Responders", style={"color": "#cc0000"},),
                                    html.P(
                                        [
                                            "Cohort: ",
                                            html.Strong("Melanoma & miraclib treatment & PBMC samples."),
                                            " Each box shows the distribution of relative frequency (%) across all"
                                            " samples in that group. Individual data points are overlaid."
                                            "Statistical significance is assessed with a two-sided ",
                                            html.Strong("Welch’s t-test"),
                                            " (α = 0.05).",
                                        ]
                                    ),
                                ],
                                style={**CARD, "marginBottom": "0", "borderBottomLeftRadius": "0", "borderBottomRightRadius": "0"},
                            ),

                            html.Div(
                                dcc.Graph(
                                    figure=make_boxplot(part3_plot_df),
                                    config={"displayModeBar": True},
                                    style={"height": "520px"},
                                ),
                                style={**CARD, "borderTopLeftRadius": "0", "borderTopRightRadius": "0", "borderTop": "1px solid #ffffff"},
                            ),

                            html.Div(
                                [
                                    html.H4("Statistical Test Results", style={"marginTop": 0}),
                                    html.P(
                                        "The Welch’s t-test compares group means without assuming equal variance, making it appropriate for comparing immune cell frequency distributions between responders and non-responders.",
                                        style={"color": "#000000", "fontSize": "13px"},
                                    ),
                                    dash_table.DataTable(
                                        data=part3_stat_df.to_dict("records"),
                                        columns=[{"name": c, "id": c} for c in part3_stat_df.columns],
                                        style_header=TABLE_HEADER,
                                        style_cell={**TABLE_CELL, "textAlign": "center"},
                                        style_data_conditional=[
                                            *TABLE_ODD,
                                            {
                                                "if": {
                                                    "filterquery": '{Significant (p < 0.05)} = "Yes"',
                                                    "column_id": "Significant (p < 0.05)",
                                                },
                                                "backgroundColor": "#ffffff",
                                                "color": "#00cc00",
                                                "fontWeight": "bold",
                                            },
                                            {
                                                "if": {
                                                    "filterquery": '{Significant (p < 0.05)} = "No"',
                                                    "column_id": "Significant (p < 0.05)",
                                                },
                                                "color": "#cc0000",
                                            },
                                        ],
                                    ),
                                ],
                                style=CARD,
                            ),
                        ],
                    ),

                    # Part 4
                    dcc.Tab(
                        label="Part 4: Baseline Subset Analysis",
                        children=[
                            html.Div(
                                [
                                    html.H3("Baseline Subset Analysis", style={"color": "#cc0000"},),
                                    html.Div(
                                        [
                                            html.P(
                                                [
                                                    "Query filter: ",
                                                    html.Code(
                                                        "condition = 'melanoma'  AND  treatment = 'miraclib'  "
                                                        "AND  sample_type = 'PBMC'  AND  time_from_treatment_start = 0",
                                                        style={"backgroundColor": "#ffffff", "padding": "2px 6px", "borderRadius": "4px"},
                                                    ),
                                                ],
                                                style={"marginBottom": "8px"},
                                            ),
                                            html.P(
                                                [
                                                    html.Strong(f"{n_samples_p4}"),
                                                    f" samples / ",
                                                    html.Strong(f"{n_subjects_p4}"),
                                                    " unique subjects",
                                                ],
                                                style={"color": "#ff0000", "fontSize": "15px", "marginBottom": 0},
                                            ),
                                        ],
                                        style={
                                            "backgroundColor": "#ffffff",
                                            "padding": "14px 18px",
                                            "borderRadius": "6px",
                                            "marginBottom": "24px",
                                        },
                                    ),

                                    # Summary tables
                                    html.Div(
                                        [
                                            # Samples per project
                                            html.Div(
                                                [
                                                    html.H4("Samples per Project", style={"marginTop": 0}),
                                                    html.P(
                                                        "Number of baseline PBMC samples from each project.",
                                                        style={"color": "#000000", "fontSize": "13px"},
                                                    ),
                                                    mini_table(proj_df, "proj-table"),
                                                ],
                                                style={**CARD, "flex": 1, "marginRight": "16px"},
                                            ),

                                            # Subjects by response
                                            html.Div(
                                                [
                                                    html.H4("Subjects by Treatment Response", style={"marginTop": 0}),
                                                    html.P(
                                                        "Unique subjects (patients) classified by whether they "
                                                        "responded to miraclib.",
                                                        style={"color": "#000000", "fontSize": "13px"},
                                                    ),
                                                    mini_table(resp_df, "resp-table"),
                                                ],
                                                style={**CARD, "flex": 1, "marginRight": "16px"},
                                            ),

                                            # Subjects by sex
                                            html.Div(
                                                [
                                                    html.H4("Subjects by Sex", style={"marginTop": 0}),
                                                    html.P(
                                                        "Unique subjects split by biological sex.",
                                                        style={"color": "#000000", "fontSize": "13px"},
                                                    ),
                                                    mini_table(sex_df, "sex-table"),
                                                ],
                                                style={**CARD, "flex": 1},
                                            ),
                                        ],
                                        style={"display": "flex", "alignItems": "flex-start", "flexWrap": "wrap"},
                                    ),
                                ],
                                style=CARD,
                            )
                        ],
                    ),
                ],
                colors={"border": "#ffffff", "primary": "#ff0000", "background": "#ffffff"},
                style={"fontSize": "14px", "fontFamily": "'Segoe UI', sans-serif"},
            ),
            style={"padding": "28px 32px", "backgroundColor": "#000000", "minHeight": "100vh"},
        ),
    ],
    style={"fontFamily": "'Segoe UI', sans-serif"},
)

@callback(
    Output("freq-table", "data"),
    Output("freq-info", "children"),
    Input("f-project", "value"),
    Input("f-condition", "value"),
    Input("f-treatment", "value"),
    Input("f-sample-type", "value"),
    Input("f-population", "value"),
)
def update_freq_table(project, condition, treatment, sample_type, population):
    df = freq_df.copy()

    if project != "__all__":
        df = df[df["project"] == project]
    if condition != "__all__":
        df = df[df["condition"] == condition]
    if treatment != "__all__":
        df = df[df["treatment"] == treatment]
    if sample_type != "__all__":
        df = df[df["sample_type"] == sample_type]
    if population != "__all__":
        df = df[df["population"] == population]

    display = df[["sample", "population", "total_count", "count", "percentage"]].copy()
    info = f"Showing {len(display):,} rows ({display['sample'].nunique():,} samples)"
    return display.to_dict("records"), info

if __name__ == "__main__":
    app.run(debug=False, port=8050)
