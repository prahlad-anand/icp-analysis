# Immune Cell Population Analysis

Interactive dashboard for exploring immune cell population data from a clinical trial.

---

## Quick Start

Run the following commands in sequence:

```
conda create --name teiko-test python=3.13                                
git clone https://github.com/prahlad-anand/icp-analysis.git
pip install -r requirements.txt
python load_data.py
python app.py
```

A link or pop-up to the dashboard will be displayed.

---

## Database Schema

### Entity-Relationship Overview

```
subjects ──< samples
```

### Table Definitions

**`subjects`** — one row per patient/subject

| Column     | Type    | Notes                                         |
|------------|---------|-----------------------------------------------|
| subject_id | TEXT    | Primary key                                   |
| project_id | TEXT    |                                               |
| condition  | TEXT    | Disease indication (melanoma, carcinoma, …)   |
| age        | INTEGER |                                               |
| sex        | TEXT    | M / F                                         |
| treatment  | TEXT    | Drug administered (miraclib, phauximab, none) |
| response   | TEXT    | yes / no / NULL (healthy controls)            |

**`samples`** — one row per biological sample, including cell counts

| Column                    | Type    | Notes                               |
|---------------------------|---------|-------------------------------------|
| sample_id                 | TEXT    | Primary key                         |
| subject_id                | TEXT    | FK → subjects                       |
| sample_type               | TEXT    | PBMC / WB                           |
| time_from_treatment_start | INTEGER | Days since treatment began (0/7/14) |
| b_cell                    | INTEGER |                                     |
| cd8_t_cell                | INTEGER |                                     |
| cd4_t_cell                | INTEGER |                                     |
| nk_cell                   | INTEGER |                                     |
| monocyte                  | INTEGER |                                     |

### Design Rationale

**Normalisation eliminates redundancy.** In the raw CSV every row repeats the subject's age, sex, condition, treatment and response for each sample. Splitting this into `subjects` and `samples` both stores each fact exactly once, and removes the risk of contradictory subject metadata across rows.

**Scalability**

| Issue | Handling |
|---|---|
| Scaling for samples         | Linear |
| Scaling for subjects        | Linear |
| Analytics at scale          | The shape of the database (multiple samples matching to a single subject) translates well across data warehousing applications such as Snowflake, should SQLite becomes a bottleneck. |
| Analytics by multiple users | PostgreSQL migration requires only a driver swap. |

---

## Code Structure

### `load_data.py`

Reads `cell-count.csv` once, deduplicates subjects using Python sets and writes to the tables in dependency order. It is idempotent too, dropping and recreating the database on each run.

### `app.py`

Three primary sections:

1. **Data-loading functions** (`build_frequency_df`, `build_stat_data`,
   `build_subset_data`) — each returns a `pandas.DataFrame` by issuing a focused SQL
   query and doing any reshaping (melt, percentage calculation, statistical tests) in
   Python. Keeping SQL simple and doing aggregation in pandas means the code is readable
   and unit-testable without a database.

2. **Layout** - a three-tab Dash layout, one tab per analytical part. Styling constants
   (`CARD`, `TABLE_HEADER`, etc.) are defined once at module level so the look can be
   changed globally.

3. **Callbacks** - Part 2 has a reactive callback (the five dropdown filters). Parts 3 and 4 produce static figures and tables that are computed once at startup, keeping the app responsive.

**Dash?** Dash produces a single-page application with server-side Python callbacks: all statistical computations stay in Python/scipy. Also renders Plotly figures natively, which allows for interactive zoom/pan/hover on the boxplot.

**Precomputation?** The melanoma & miraclib & PBMC cohort filter and the baseline subset filter are fixed. Running statistical tests and subset aggregations once at startup avoids recomputation on page load.

---

## Dashboard

The dashboard is served locally at:

```
http://127.0.0.1:8050
```
