
with open("clean_genes.tsv","w") as output:
    with open('annot.gff','r') as read_in:
        for line in read_in:
            l = line.strip('\n').split('\t')
            if "#" not in line:
                if l[2] == "gene":
                    start = str(l[3])
                    end = str(l[4])
                    orientacion = str(l[6])
                    info = l[8].split(";")[0].replace("ID=gene-","")
                    sentence = info+"\t"+orientacion+"\t"+start+"\t"+end+"\n"
                    output.write(sentence)

