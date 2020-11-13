import cobra.test
model=cobra.io.load_json_model("/home/agustin/FBA_Tesis/Paper_IEK1011/supplementaryMaterial/12918_2018_557_MOESM3_ESM/iEK1011_m7H10_media.json")
cobra.io.save_matlab_model(model,"iEK1011_m7H10_media.mat")
