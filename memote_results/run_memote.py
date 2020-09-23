import cobra
import memote

model=cobra.io.read_sbml_model("iEK1008.xml")
result=memote.test_model(model, skip=["test_stoichiometric_consistency"], results=True)
report = memote.snapshot_report(result[1], config=None, html=True)
with open("report.html", "w") as handle:
    handle.write(report)
