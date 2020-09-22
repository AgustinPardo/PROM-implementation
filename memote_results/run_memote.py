import cobra
import memote

model=cobra.io.read_sbml_model("iEK1008.xml")
result=memote.test_model(model, skip=["test_consistency"], results=True)
memote.snapshot_report(result[1], config=None, html=True)
