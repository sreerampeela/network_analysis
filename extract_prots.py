from Bio import SeqIO

fastas1 = SeqIO.parse("373153_prots.fa","fasta")
fastas2 = SeqIO.parse("171101_prots.fa","fasta")
fastas3 = SeqIO.parse("170187_prots.fa","fasta")

filein = open("ncbi_3_blast.tab",'r',newline="\n",encoding="utf-8")

prot_list1 = [i.split("\t")[0].rstrip() for i in filein.readlines()]
prot_list2 = [i.split("\t")[1].rstrip() for i in filein.readlines()]
prot_list = list(set(prot_list1+prot_list2))
print(len(prot_list))
print(prot_list)

records = []
# print(len(prot_list1))
for rec1 in fastas1:
    if rec1.id in prot_list:
        # print(prot_list.index(rec.id))
        records.append(rec1)

for rec2 in fastas2:
    if rec2.id in prot_list:
        # print(prot_list.index(rec.id))
        records.append(rec2)

for rec3 in fastas3:
    if rec3.id in prot_list:
        # print(prot_list.index(rec.id))
        records.append(rec3)

print(len(records))

SeqIO.write(records,"pan_proteome_spn.fasta","fasta")

filein.close()