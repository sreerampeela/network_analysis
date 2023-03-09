import pandas as pd
blast_res = open("spn_proteome_blast.csv",'r',newline='\n',encoding='utf-8')

essential_pan = []
for i in list(blast_res.readlines()):
    # print(i)
    if i.split(",")[6] == "yes" and i.split(",")[8] == "yes" and i.split(",")[10] == "yes":
        # print(i.split(",")[0])
        essential_pan.append(i.split(",")[0])

print(len(essential_pan))


blast_res.close()

df = pd.read_csv("ncbi_3_blast.tab",sep="\t",header=None)


# df["source_essential"] = [None] * len(df[0])
# df["target_essential"] = [None] * len(df[0])

print(df.head())
sources = df[0]
targets = df[1]
source_essential = []
target_essential = []
for i in sources:
    if i in essential_pan:
        source_essential.append("y")
    else:
        source_essential.append("n")

for j in targets:
    if j in essential_pan:
        target_essential.append("y")
    else:
        target_essential.append("n")

df["source_essential"] = source_essential
df["target_essential"] = target_essential

edge_type = []
for i in range(len(source_essential)):
    if source_essential[i] == "n" and target_essential[i] == "n":
        edge_type.append("N-N")
    elif source_essential[i] == "n" and target_essential[i] == "y":
        edge_type.append("N-E")
    elif source_essential[i] == "y" and target_essential[i] == "n":
        edge_type.append("E-N")
    elif source_essential[i] == "y" and target_essential[i] == "y":
        edge_type.append("E-E")

# print(edge_type)

df["edge_type"] = edge_type

print(df)
df.to_csv("essential_appended.csv")

# source	target	gene_source	gene_target	taxon	score	neighborhood	fusion	phylogen_profile	co_expression	expt	db_score	text_mining	source_essential	target_essential

