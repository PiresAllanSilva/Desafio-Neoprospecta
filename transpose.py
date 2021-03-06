import pandas as pd
data = pd.read_csv("C:/Users/karen/OneDrive/Área de Trabalho/Állan/Neoprospecta/tables/otu_table_tax_amostras.tsv", sep ="\\t", header = 0)
#time = pd.read_csv("C:/Users/karen/OneDrive/Área de Trabalho/Állan/Neoprospecta/metadata.csv", sep ="\\t")

#print(time)
df = pd.DataFrame(data, columns = ["F3D0_S188","F3D141_S207","F3D142_S208","F3D143_S209","F3D144_S210","F3D145_S211",
                                   "F3D146_S212","F3D147_S213","F3D148_S214","F3D149_S215","F3D150_S216","F3D1_S189",
                                   "F3D2_S190","F3D3_S191","F3D5_S193","F3D6_S194","F3D7_S195","F3D8_S196","F3D9_S197"])
df = df.transpose()
df.to_csv("C:/Users/karen/OneDrive/Área de Trabalho/Állan/Neoprospecta/otu.csv", index = True)