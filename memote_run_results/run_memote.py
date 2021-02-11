import cobra
import memote

model=cobra.io.read_sbml_model("/home/agustin/FBA_Tesis/H37Rv-CHROM.sbml")
result=memote.test_model(model, skip=["test_consistency"], results=True)
report = memote.snapshot_report(result[1], config=None, html=True)
with open("H37Rv-CHROM.html", "w") as handle:
    handle.write(report)
