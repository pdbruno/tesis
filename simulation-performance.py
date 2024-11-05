import pandas as pd


def latex_scientific_notation(f):
    [base, exp] = "{:.2E}".format(f).split("E")
    return "$" + base + "\\times10^{" + exp + "}$"


index = ["chA", "chB", "d", "pA", "pB"]
simulacion = pd.read_csv("haar-avg-sim.csv").set_index(index).sort_index()
teorico = pd.read_csv("teorico.csv").set_index(index).sort_index()
teorico["score"] = (
    teorico["score"]
    .str.replace("1/2", "0.5")
    .str.replace("2/3", "0.6666666666666666")
    .str.replace("1/Sqrt[2]", "0.7071067811865475")
    .str.replace("Pi/4", "0.7853981633974483")
    .values.astype("float64")
)
simulacion["score_teorico"] = teorico["score"].values
grouped_diffs = (
    (simulacion["score"] - simulacion["score_teorico"])
    .abs()
    .groupby(["chA", "chB", "d"])
)
final = grouped_diffs.agg(["mean", "std"]).reset_index()
final.to_latex(
    "simulation-performance.txt",
    float_format=latex_scientific_notation,
    formatters={},
    index=False,
)
print(final[final["mean"].abs() == final["mean"].abs().min()])
print(final[final["mean"].abs() == final["mean"].abs().max()])