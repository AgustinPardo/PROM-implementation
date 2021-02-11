import cobra
import memote

model1=cobra.io.load_json_model("iNJ661.json")
model2=cobra.io.load_json_model("/home/agustin/FBA_Tesis/Paper_IEK1011/supplementaryMaterial/12918_2018_557_MOESM3_ESM/iEK1011_m7H10_media.json")
result=memote.suite.reporting.diff(model1, model2, skip=["test_consistency"], results=True)
report = memote.snapshot_report(result[1], config=None, html=True)
with open("iNJ661.html", "w") as handle:
    handle.write(report)
