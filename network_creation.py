import os
import networkx as nx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests


def read_network(infile, keep_scores=True, cutoff=700):
    """To read a network file from STRINGdb and convert into a graph object"""
    with open(infile, 'r', encoding='utf-8', newline='\n') as fin:
        data = [i.rstrip() for i in fin.readlines()]
    edge_data = []
    for i in data[1:]:
        source_node = i.split(" ")[0].split(".")[1]
        # print(source_node)
        target_node = i.split(" ")[1].split(".")[1]
        score = i.split(" ")[-1]
        if float(score) >= cutoff:
            # print(type(score))
            if keep_scores is True:
                p = "\t".join([source_node, target_node, score])
            else:
                p = "\t".join([source_node, target_node])
            edge_data.append(p)
    if keep_scores is True:
        graph_name = nx.parse_edgelist(
            edge_data, delimiter="\t", create_using=nx.DiGraph, data=(("weight", float),))
    else:
        graph_name = nx.parse_edgelist(
            edge_data, delimiter="\t", create_using=nx.DiGraph)

    return graph_name


def run_blast(db_file, file1, file2):
    """Runs blastp on a query proteome and a reference proteome 'db_file'. Only one best hit per query sequence returned. 
    The results are exported to a csv file 'file2'."""

    # creating a blast database
    mkdb_cmd = f"makeblastdb -in {db_file} -dbtype prot -title {db_file} -logfile blastp_run.log"
    os.system(mkdb_cmd)

    # running blastp with the database built above
    cmd = f"blastp -db {db_file} -query {file1}" \
        f" -outfmt '6 delim=, qseqid sseqid slen qlen length pident evalue score qstart qend sstart send gaps positive' " \
        f"-out {file2} -num_alignments 1 -num_threads 4 -mt_mode 1 "
    os.system(cmd)


def create_ref_proteome(file1, file2, blast_res):
    """Merge two proteomes basesd on sequence similarity between the two. File 2 is considered as query proteome. 
    Uses results (csv file) from blastp"""
    # column_names = ["qseqid", "sseqid", "slen", "qlen", "length", "pident",
    #                "evalue", "score", "qstart", "qend", "sstart", "send", "gaps", "positive"]
    infile1 = open(file1, 'r', encoding='utf-8')
    infile2 = open(file2, 'r', encoding='utf-8')
    file_1 = SeqIO.parse(infile1, "fasta")
    file_2 = SeqIO.parse(infile2, "fasta")

    temp_prot = open("temp_proteome4.fasta", 'a')

    records_1 = [i for i in file_1]
    records_2 = [j for j in file_2]
    print(records_1)

    blast_results = open(blast_res, 'r', encoding='utf-8', newline='\n')
    blast_data = blast_results.readlines()

    for i in blast_data:
        query_seq = i.split(",")[0].split(".")[1]
        best_hit = i.split(",")[1].split(".")[1]
        percent_iden = float(i.split(",")[5])
        subj_cov = round(float(i.split(",")[4])/float(i.split(",")[2])*100, 2)
        query_cov = round(float(i.split(",")[4])/float(i.split(",")[3])*100, 2)
        refseq = file1.split("_")[0] + "." + best_hit
        # print(refseq)
        alt_seq = file2.split("_")[0] + "." + query_seq
        if (percent_iden >= 95) and (subj_cov >= 90) and (query_cov >= 90):
            print(
                f"match found between {query_seq} and {best_hit} with {percent_iden} identity and {query_cov} query coverage and {subj_cov} subject coverage")
            for rec1 in records_1:
                if rec1.id == refseq:
                    sequence = rec1.seq
                    seq_ref = SeqRecord(id=rec1.id, seq=sequence)
                    SeqIO.write(seq_ref, temp_prot, "fasta")
                    # print(seq_ref)

        else:
            seq1 = refseq
            seq2 = alt_seq
            for rec1 in records_1:
                if rec1.id == seq1:
                    sequence = rec1.seq
                    seq_ref = SeqRecord(id=rec1.id, seq=sequence)
                    SeqIO.write(seq_ref, temp_prot, "fasta")
                    print(seq_ref)
            for rec2 in records_2:
                if rec2.id == seq2:
                    sequence2 = rec2.seq
                    seq_alt = SeqRecord(id=rec2.id, seq=sequence2)
                    SeqIO.write(seq_alt, temp_prot, "fasta")
                    print(seq_alt)

    blast_results.close()
    temp_prot.close()


def get_pangeneids(blast_res):
    """Get gene ids from BLAST results based on similarity and coverage thresholds"""
    blast_results = open(blast_res, 'r', encoding='utf-8', newline='\n')
    gene_ids = []
    for blast_data in blast_results.readlines():
        query_seq = blast_data.split(",")[0].split(".")[1]
        best_hit = blast_data.split(",")[1].split(".")[1]
        percent_iden = float(blast_data.split(",")[5])
        subj_cov = round(float(blast_data.split(
            ",")[4])/float(blast_data.split(",")[2])*100, 2)
        query_cov = round(float(blast_data.split(
            ",")[4])/float(blast_data.split(",")[3])*100, 2)
        if (percent_iden >= 95) and (subj_cov >= 90) and (query_cov >= 90):
            gene_ids.append(best_hit)
        else:
            gene_ids.append(best_hit)
            gene_ids.append(query_seq)
    return gene_ids


def get_string_network(idslist, foutput):
    """Imports the STRING network via API for the sequences in fname File"""
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"
    request_url = "/".join([string_api_url, output_format, method])

    # records = SeqIO.parse(fname, "fasta")
    for record in idslist:
        string_id = record
        print(string_id)
        params = {
            "limit": 20,
            "identifiers": string_id,
            "species": 2,
            "caller_identity": "srirambds.murthy@gmail.com",
            "required_score": 700
        }
        response = requests.post(request_url, data=params)
        with open(foutput, 'a') as fout:
            fout.write(response.text)
        print(f"Written protein interactions for {string_id}")
        os.system("sleep 3")


def get_BLAST_hits_dict(blast_res1, blast_res2):
    """Input BLAST results file. Creates a dictionary with reference gene id as key and orthologs as values"""
    ref_ids_keys = []
    result_dict = dict()
    with open(blast_res1, 'r', encoding='utf-8') as fin1:
        for i in fin1.readlines():
            refid = i.split(",")[1].split(".")[1]
            altid = i.split(",")[0].split(".")[1]
            percent_iden = float(i.split(",")[5])
            subj_cov = round(
                float(i.split(",")[4])/float(i.split(",")[2])*100, 2)
            query_cov = round(
                float(i.split(",")[4])/float(i.split(",")[3])*100, 2)
            if (percent_iden >= 95) and (subj_cov >= 90) and (query_cov >= 90):
                ref_ids_keys.append(refid)
                datadict = {altid: refid}
                # print(datadict)
                result_dict.update(datadict)

    with open(blast_res2, 'r', encoding='utf-8') as fin2:
        for j in fin2.readlines():
            refid2 = j.split(",")[1].split(".")[1]
            altid2 = j.split(",")[0].split(".")[1]
            percent_iden = float(j.split(",")[5])
            subj_cov = round(
                float(j.split(",")[4])/float(j.split(",")[2])*100, 2)
            query_cov = round(
                float(j.split(",")[4])/float(j.split(",")[3])*100, 2)
            if (percent_iden >= 95) and (subj_cov >= 90) and (query_cov >= 90):
                datadict2 = {altid2: refid2}
                result_dict.update(datadict2)
    return result_dict


def renameNodes(file1, hits_dict):
    """Reads a network file and rename the nodes to reference gene ids only if gene ids start with SP_ or spr_. 
    Uses the blast hits dictionary created before. Returns an edgelist to be written to a file"""
    with open(file1, 'r', encoding="utf-8") as fin:
        edges = [(i.split("\t")[2], i.split("\t")[3]) for i in fin.readlines()]
# # print(len(edges))
    new_edges = edges[:]
    for edge in edges[201:]:
        source, target = edge
        if (source.startswith("SP_") is True) or (source.startswith("spr_") is True):
            if source in list(hits_dict.keys()):
                new_source = hits_dict[source]
                print(new_source, target)
                new_edges.remove(edge)
                new_edges.append((new_source, target))
    new_edges2 = set(new_edges)
    return new_edges2
