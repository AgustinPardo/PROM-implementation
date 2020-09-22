import cobra

model_Bigg = cobra.io.load_matlab_model("/home/agustin/FBA_Tesis/trabajo/iEK1008.mat")
model_m7H10 = cobra.io.load_json_model("/home/agustin/FBA_Tesis/Paper_IEK1011/supplementaryMaterial/12918_2018_557_MOESM3_ESM/iEK1011_m7H10_media.json")
model_inVivo= cobra.io.load_json_model("/home/agustin/FBA_Tesis/Paper_IEK1011/supplementaryMaterial/12918_2018_557_MOESM3_ESM/iEK1011_inVivo_media.json")
model_Bigg.name="Bigg"
model_m7H10.name="m7H10"
model_inVivo.name="inVivo"

def FBA(model, biomass_name):
    print("Model name: "+model.name)
    print("Medio:")
    print(model.medium)

    print("Biomass composition reaction:")

    print(model.reactions.get_by_id(biomass_name).reaction)

    solution=model.optimize()

    print("FBA summary:")
    print(model.summary())
    #print(solution.objective_value)

    return 0

FBA(model_Bigg,"BIOMASS__2")
print("*"*20)
FBA(model_m7H10,"biomass")
print("*"*20)
FBA(model_inVivo,"biomass")