from statsmodels.stats.proportion import proportion_confint
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import date
import datetime
import csv
import baltic as bt
import collections
# from Bio.Alphabet import IUPAC
# from pySankey.sankey import sankey
import math
import collections
import datetime as dt

import os
import csv
from Bio import SeqIO
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
import collections
import numpy as np
import pandas as pd
import datetime
from datetime import date
# import skbio
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.lines as mlines
import collections
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from matplotlib.collections import PatchCollection
from epiweeks import Week, Year
font = {'family' : 'Helvetica',
        'weight' : 'bold',
        'size'   : 18}
from matplotlib.lines import Line2D

mpl.rcParams.update({'font.size': 18})

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
mpl.rcParams.update(new_rc_params)

plt.rcParams['font.family'] = 'Helvetica'


def get_node_states_all_sites(directory,state_file,alignment):
    
    #returns a dict keys off 1-based positions in the genome
    #and the value is a list of tuples (id, base) at a given node
    # allows you to look up for a given site what the base is for a
    #given internal node or tip
    
    node_states = collections.defaultdict(list)
    c = 0
    
    ## first the reconstructed nodes
    with open(f"{directory}/{state_file}","r") as f:
        for l in f:

            if not l.startswith("#"):
                c+=1
                try:
                    node,site,state,probA,probC,probG,probT = l.rstrip("\n").split("\t")
                except:
                    print(l)
                    break
                if node != "Node":
                    if state not in ["N","-"]:
                        node_states[site].append((node,state))
                    else:
                        node_states[site].append((node,""))
    ## now the tips
    for record in SeqIO.parse(f"{directory}/{alignment}","fasta"):
        for site in node_states:
            index = int(site)-1
            base = record.seq[index]
            if base in ["T","C","A","G"]:
                node_states[site].append((record.id,base))
            else:
                node_states[site].append((record.id,""))
                
    return node_states

def get_header_str(dict_values):
    header_str = ""
    for i in sorted(dict_values, key = lambda i : i[0]):
        header_str += f"{i[0]},"
    header_str = header_str.rstrip(",")
    return header_str
    
    
def find_what_sites_vary_unambiguously(node_states,outfile):
    header_str = get_header_str(node_states["1"])
    
    with open(outfile,"w") as fw:
        fw.write(f"site,{header_str}\n")

        for site in node_states:
            info = node_states[site]
            
            # get the set of unique bases at a given site
            count = set([i[1] for i in info if i[1]])
            
            #if there's more than one
            if len(count)>1:
                
                #needs to be kep consistent with header str
                info = sorted(info, key = lambda i : i[0])
                base_str = ""
                for i in info:
                    base_str += f"{i[1]},"
                    
                base_str = base_str.rstrip(",")
                fw.write(f"{site},{base_str}\n")
    
def load_unambiguous_varying_sites(infile):
    node_states_diff = collections.defaultdict(dict)
    with open(infile,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site = row["site"]
            for col in row:
                if col != "site":
                    node_states_diff[row["site"]][col] = row[col]
    return node_states_diff

def map_site_changes_to_branches(treefile, outfile,node_states,node_states_diff): 
    my_tree=bt.loadNewick(treefile,absoluteTime=False)
    last_node = ""
    current_node = ""

    with open(outfile,"w") as fw:
        fw.write("parent,child,site,snp,dimer\n")

        for k in my_tree.Objects:
            if k.branchType == 'leaf':
                current_node = k
                current_node.traits["label"]=k.name
            else:
                current_node = k

            if last_node:
                node_name = current_node.traits["label"]
                parent_name = current_node.parent.traits["label"]
                snps = []
                for site in node_states_diff:
                    node_base = node_states_diff[site][node_name]
                    parent_base = node_states_diff[site][parent_name]

                    if node_base != parent_base:
                        if node_base in ["A","C","G","T"] and parent_base in ["A","C","G","T"]:
                            snp = f"{parent_base}->{node_base}"
                            snps.append(snp)
                            if snp == "G->A":
                                dimer_site = f"{int(site)+1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{parent_base}{dimer_base}"
                            elif snp == "C->T":
                                dimer_site = f"{int(site)-1}"
                                dimer_base = ""

                                for i in node_states[dimer_site]:
                                    if i[0] == parent_name:
                                        dimer_base = i[1]
                                dimer = f"{dimer_base}{parent_base}"
                            else:
                                dimer = ""
                            fw.write(f"{parent_name},{node_name},{site},{snp},{dimer}\n")

            last_node = current_node

def read_in_branch_snps(branch_snps):
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            branch = f"{row['parent']}_{row['child']}"
            branch_snps_dict[branch].append((row['site'],row['snp'],row['dimer'])) 
    return branch_snps_dict

def get_branch_snps_sites(branch_snps):
    all_snps = collections.Counter()
    branch_snps_dict = collections.defaultdict(list)
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            all_snps[int(row["site"])]+=1

            branch_snps_dict[int(row['site'])].append([row['parent'],row['child'],row['snp'],row['dimer']])
    
    homoplasies = {}
    for k in all_snps:
        if all_snps[k]>1:
            homoplasies[k] = all_snps[k]
            
    print(len(homoplasies))
    print(homoplasies)
    return branch_snps_dict,homoplasies
    

def get_acc_to_metadata_map(metadata):
    acc_dict = {}
    with open(metadata,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                acc_dict[row["accession"]] = row
            except:
                try:
                    acc_dict[row["\ufeffaccession"]] = row
                except:
                    pass
                
    return acc_dict

def make_reconstruction_tree_figure(outfile,branch_snps,treefile,width,height):
    
    branch_snps_dict = read_in_branch_snps(branch_snps)
    print(len(branch_snps_dict),"snps")
#     acc_dict = get_acc_to_metadata_map(metadata)
#     print(len(acc_dict))        
    
    my_tree=bt.loadNewick(treefile,absoluteTime=False)

    fig,ax = plt.subplots(figsize=(width,height),facecolor='w')

    x_attr=lambda k: k.height ## x coordinate of branches will be absoluteTime attribute
    su_func=lambda k: 50-30*k.height/my_tree.treeHeight ## size of tips
    s_func=lambda k: 50-20*k.height/my_tree.treeHeight ## size of tips
    c_func=lambda k: "dimgrey"

    increment = my_tree.treeHeight/150
    my_tree.plotTree(ax,x_attr=x_attr) ## plot branches
    my_tree.plotPoints(ax,size=s_func,colour=c_func,x_attr=x_attr) ## plot circles at tips
    mpl.rcParams['font.family'] = 'sans-serif'
    

    
    
    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits["label"]=k.name

        node_name = current_node.traits["label"]
        try:
            parent_name = current_node.parent.traits["label"]
        except:
            continue
        branch_name= f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
#                 print(branch_name, len(branch_snps_dict[branch_name]))
            snp_placement = current_node.parent.height + increment
            for s in branch_snps_dict[branch_name]:
                site,snp,dimer = s
                if snp == "G->A":
                    if dimer in ["GA"]:
                        snps.append((1,"#995e62"))
                    else:
                        snps.append((2,"#d9b660"))
                elif snp == "C->T":
                    if dimer in ["TC"]:
                        snps.append((1,"#995e62"))
                    else:
                        snps.append((2,"#d9b660"))
                else:
                    snps.append((2,"#d9b660"))

            for snp in sorted(snps, key = lambda x : x[0]):
                plt.scatter([snp_placement],[k.y+0.5],color=snp[1],s=30)
                snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    ax.tick_params(axis='y',size=0)
    ax.tick_params(axis='x',size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.savefig(f"figures/{outfile}.svg");
    plt.savefig(f"figures/{outfile}.png",bbox_inches='tight', 
                   transparent=True);
    plt.show()
    

def make_reconstruction_tree_figure_w_labels(outfile,branch_snps,treefile,width,height):
    
    branch_snps_dict = read_in_branch_snps(branch_snps)
    print(len(branch_snps_dict),"branches")
#     acc_dict = get_acc_to_metadata_map(metadata)
#     print(len(acc_dict))        
    
    my_tree=bt.loadNewick(treefile,absoluteTime=False)

    fig,ax = plt.subplots(figsize=(width,height),facecolor='w')

    x_attr=lambda k: k.height ## x coordinate of branches will be absoluteTime attribute
    su_func=lambda k: 50-30*k.height/my_tree.treeHeight ## size of tips
    s_func=lambda k: 50-20*k.height/my_tree.treeHeight ## size of tips
    c_func=lambda k: "dimgrey"

    increment = my_tree.treeHeight/150
    my_tree.plotTree(ax,x_attr=x_attr) ## plot branches
    my_tree.plotPoints(ax,size=s_func,colour=c_func,x_attr=x_attr) ## plot circles at tips
    mpl.rcParams['font.family'] = 'sans-serif'
    
    
    target_func=lambda k: k.is_leaf() ## which branches will be annotated
    text_func=lambda k: k.name ## what text is plotted
    text_x_attr=lambda k: k.height+(increment*4) ## where x coordinate for text is

    my_tree.addText(ax,x_attr=text_x_attr,target=target_func,text=text_func) #
    
    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits["label"]=k.name

        node_name = current_node.traits["label"]
        try:
            parent_name = current_node.parent.traits["label"]
        except:
            continue
        branch_name= f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
#                 print(branch_name, len(branch_snps_dict[branch_name]))
            snp_placement = current_node.parent.height + increment
            for s in branch_snps_dict[branch_name]:
                site,snp,dimer = s
                if snp == "G->A":
                    if dimer in ["GA"]:
                        snps.append((1,"#995E62"))
                    else:
                        snps.append((3,"#D9B660"))
                elif snp == "C->T":
                    if dimer in ["TC"]:
                        snps.append((2,"#995E62"))
                    else:
                        snps.append((3,"#D9B660"))
                else:
                    snps.append((4,"#D9B660"))

            for snp in sorted(snps, key = lambda x : x[0]):
                plt.scatter([snp_placement],[k.y+0.5],color=snp[1],s=30)
                snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top','right','left']]
    ax.tick_params(axis='y',size=0)
#     ax.tick_params(axis='x',size=0)

    ax.set_yticklabels([])
#     ax.set_xticklabels([])

    plt.savefig(f"figures/{outfile}.svg");
    plt.savefig(f"figures/{outfile}.png",bbox_inches='tight', 
                   transparent=True);
    plt.show()
    
def generate_reconstruction_files(directory, alignment,treefile):
    
    node_states = get_node_states_all_sites(directory,f"{alignment}.state",alignment)
        
    find_what_sites_vary_unambiguously(node_states,f"{directory}/{treefile}.state_differences.csv")
    
def load_info(directory, alignment,treefile,treefigureout,width=25,height=30):
    node_states = get_node_states_all_sites(directory,f"{alignment}.state",alignment)
    node_states_diff = load_unambiguous_varying_sites(f"{directory}/{treefile}.state_differences.csv")

    map_site_changes_to_branches(f"{directory}/{treefile}",
                                 f"{directory}/{treefile}.branch_snps.reconstruction.csv",
                                 node_states,
                                 node_states_diff)

    make_reconstruction_tree_figure(treefigureout,
                                    f"{directory}/{treefile}.branch_snps.reconstruction.csv",
                                    f"{directory}/{treefile}",width,height)
    
def make_apobec_context_mutation_count_plot(branch_snps,out_counts,outfigure):
    type_counter = collections.Counter()
    apobec_counter = collections.Counter()
    all_snps_count = 0
    all_apobec_count = 0
    with open(branch_snps,"r") as f:
        
        reader = csv.DictReader(f)
        for row in reader:
#             print(row)
            if row["parent"] != "Node1":
                if row["dimer"] in ["TC","GA"]:

                    apobec_counter[f'{row["snp"]}']+=1
                    all_apobec_count +=1

                type_counter[f'{row["snp"]}']+=1
                all_snps_count +=1
                
    print(all_snps_count, all_apobec_count, all_apobec_count/all_snps_count)  
    with open(out_counts,"w") as fw:
        fw.write("snp,count,target\n")            
        for i in sorted(type_counter, key=lambda x : type_counter[x], reverse=True):
            if i in apobec_counter:
                fw.write(f"{i},{type_counter[i]},{apobec_counter[i]}\n")
            else:
                fw.write(f"{i},{type_counter[i]},0\n")
    
    df = pd.read_csv(out_counts)
    fig,ax= plt.subplots(figsize=(6,4),facecolor='w',frameon=False)

    sns.barplot(x = 'snp', y = 'count', data = df, color = 'dimgrey',alpha=0.4)

    sns.barplot(x = 'snp', y = 'target', data = df, color = 'steelblue',alpha=0.6)

    [ax.spines[loc].set_visible(False) for loc in ['top','right']]
    plt.xlabel("SNP")
    plt.ylabel("Count")

    plt.savefig(f"figures/{outfigure}.svg");
    plt.savefig(f"figures/{outfigure}.png",bbox_inches='tight', 
                   transparent=True);

    plt.show();

def make_apobec_context_internal_plot(branch_snps,out_counts,outfigure):
    type_counter = collections.Counter()
    apobec_counter = collections.Counter()
    all_snps_count = 0
    all_apobec_count = 0
    with open(branch_snps,"r") as f:
        
        reader = csv.DictReader(f)
        for row in reader:
#             print(row)
            if row["child"].startswith("Node") and row["parent"] != "Node1":
                if row["dimer"] in ["TC","GA"]:

                    apobec_counter[f'{row["snp"]}']+=1
                    all_apobec_count +=1

                type_counter[f'{row["snp"]}']+=1
                all_snps_count +=1
                
    print(all_snps_count, all_apobec_count, all_apobec_count/all_snps_count)  
    with open(out_counts,"w") as fw:
        fw.write("snp,count,target\n")            
        for i in sorted(type_counter, key=lambda x : type_counter[x], reverse=True):
            if i in apobec_counter:
                snp = i.replace("->","")
                fw.write(f"{snp},{type_counter[i]},{apobec_counter[i]}\n")
            else:
                snp = i.replace("->","")
                fw.write(f"{snp},{type_counter[i]},0\n")
    
    df = pd.read_csv(out_counts)
    fig,ax= plt.subplots(figsize=(8,4),facecolor='w',frameon=False)

    sns.barplot(x = 'snp', y = 'count', data = df, color = 'dimgrey',alpha=0.4)

    sns.barplot(x = 'snp', y = 'target', data = df, color = 'steelblue',alpha=0.6)

    [ax.spines[loc].set_visible(False) for loc in ['top','right']]
    plt.xlabel("SNP")
    plt.ylabel("Count")

    plt.savefig(f"figures/{outfigure}.svg");
    plt.savefig(f"figures/{outfigure}.png",bbox_inches='tight', 
                   transparent=True);

    plt.show();
    
    
def make_apobec_context_external_plot(branch_snps,out_counts,outfigure):
    type_counter = collections.Counter()
    apobec_counter = collections.Counter()
    all_snps_count = 0
    all_apobec_count = 0
    with open(branch_snps,"r") as f:
        
        reader = csv.DictReader(f)
        for row in reader:
#             print(row)
            if not row["child"].startswith("Node") and row["parent"] != "Node1":
                if row["dimer"] in ["TC","GA"]:

                    apobec_counter[f'{row["snp"]}']+=1
                    all_apobec_count +=1

                type_counter[f'{row["snp"]}']+=1
                all_snps_count +=1
                
    print(all_snps_count, all_apobec_count, all_apobec_count/all_snps_count)  
    with open(out_counts,"w") as fw:
        fw.write("snp,count,target\n")            
        for i in sorted(type_counter, key=lambda x : type_counter[x], reverse=True):
            if i in apobec_counter:
                fw.write(f"{i},{type_counter[i]},{apobec_counter[i]}\n")
            else:
                fw.write(f"{i},{type_counter[i]},0\n")
    
    df = pd.read_csv(out_counts)
    fig,ax= plt.subplots(figsize=(12,4),facecolor='w',frameon=False)

    sns.barplot(x = 'snp', y = 'count', data = df, color = 'dimgrey',alpha=0.4)

    sns.barplot(x = 'snp', y = 'target', data = df, color = 'steelblue',alpha=0.6)

    [ax.spines[loc].set_visible(False) for loc in ['top','right']]
    plt.xlabel("SNP")
    plt.ylabel("Count")

    plt.savefig(f"figures/{outfigure}.svg");
    plt.savefig(f"figures/{outfigure}.png",bbox_inches='tight', 
                   transparent=True);

    plt.show();
    
def get_aa_position(index,based=0):
    if based != 0:
        index -= based 

    remainder = index%3
    position_dict = {1:2,2:3,0:1}
    
    return position_dict[remainder]

def reverse_aa_position(start,end,site):
    index_dict = {}
    index = 0
    rev_position_dict = {1:2,2:3,0:1}
    for i in reversed(range(start,end)):
        remainder = index%3
        index_dict[i] = rev_position_dict[remainder]
        index +=1 

    position = index_dict[site]
    
    return position

def get_codon_indexes(aa_position,index):
    codon = []
    if aa_position == 1:
        codon = [index,index+1,index+2]
    elif aa_position == 2:
        codon = [index-1,index,index+1]
    elif aa_position == 3:
        codon =  [index-2,index-1,index]
    else:
        print("incorrect aa position")
        
    return codon

def get_codon_indexes_rev_strand(position,index):
    codon = []
    if position == 3:
        codon = [index,index+1,index+2]
    elif position == 2:
        codon = [index-1,index,index+1]
    elif position == 1:
        codon =  [index-2,index-1,index]
    return codon


def get_gene_boundaries():
    genes = {}
    gene_id = 0
    with open("/Users/s1680070/repositories/alignHPXV/squirrel/data/gene_boundaries.csv","r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene_id +=1
            name = f"{row['Name'].replace(' ','_')}_{gene_id}"
            start = int(row["Minimum"]) - 1
            end = int(row["Maximum"])
            length = int(row["Length"])
            direction = row["Direction"]
            genes[(start,end)]=(name,length,direction)
    return genes

def get_grantham_scores():
    grantham_scores = {}

    with open("metadata/grantham_score.txt","r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for col in row:
                if col!="FIRST":
                    mutation = f"{row['FIRST']}{col}"

                    if row[col] != "0":
                        grantham_scores[mutation] = int(row[col])                    
                        grantham_scores[mutation[::-1]] = int(row[col])
    return grantham_scores
 
    
def categorise_amino_acid_mutation(aa1,aa2,grantham_scores):
    
    mutation_category = ""
    if aa1 == aa2:
        mutation_category = "synonymous"
    else:
        if aa2 == '*':
            mutation_category = "nonsense"
        else:
            mutation_category = "nonsynonymous"

    if f"{aa1}{aa2}" in grantham_scores:
        score = grantham_scores[f"{aa1}{aa2}"]
        if score < 51:
            prediction = "conservative"
        elif score <101:
            prediction = "moderately conservative"
        elif score <151:
            prediction = "moderately radical"
        else:
            prediction = "radical"
    else:
        score = "NA"
        prediction = "NA"
    
    return mutation_category,score,prediction

def reconstruct_amino_acid_mutations(branch_snps,node_states,outfile):
    branch_snps_dict,homoplasies = get_branch_snps_sites(branch_snps)
    genes = get_gene_boundaries()
    grantham_scores = get_grantham_scores()
    
    fw = open(outfile,"w")
    fw.write("site,gene,direction,snp,dimer,apobec,aa_position,parent,parent_codon,parent_aa,")
    fw.write("child,child_codon,child_aa,mutation_category,score,prediction,homoplasy,occurrence\n")
    
    for site in branch_snps_dict:
        
        homoplasy = "False"
        occurrence = "1"
        if site in homoplasies:
            homoplasy = "True"
            occurrence = f"{homoplasies[site]}"
        site_found = False
        for gene in genes:
            start,end=gene

            if site in range(start,end):
                
                site_found = True
                name,length,direction = genes[gene]
                for site_snp in branch_snps_dict[site]:
                    parent,child,snp,dimer = site_snp

                    if direction == "forward":
                        aa_position = get_aa_position(site,gene[0])
                        codon_indexes = get_codon_indexes(aa_position,site)
                    else:
                        aa_position = reverse_aa_position(start,end,site)
                        codon_indexes = get_codon_indexes_rev_strand(aa_position,site)

                    parent_codon = []
                    child_codon = []

                    for base in codon_indexes:
                        reconstruction = node_states[f"{base}"]
                        for node in reconstruction:
                            if node[0] == parent:
                                parent_codon.append(node[1])

                            elif node[0] == child:
                                child_codon.append(node[1])

                    parent_codon = "".join(parent_codon)
                    child_codon = "".join(child_codon)

                    parent_codon = Seq(parent_codon)
                    child_codon = Seq(child_codon)

                    if direction == "reverse":
                        parent_codon= parent_codon.reverse_complement()
                        child_codon= child_codon.reverse_complement()

                    parent_aa = parent_codon.translate()
                    child_aa = child_codon.translate()

                    mutation_category,score,prediction = categorise_amino_acid_mutation(parent_aa,child_aa,grantham_scores)
                    apobec = "False"
                    if snp in ["C->T","G->A"] and dimer in ["GA","TC"]:
                        apobec = "True"

                    fw.write(f"{site},{name},{direction},{snp},{dimer},{apobec},{aa_position},{parent},{parent_codon},{parent_aa},{child},{child_codon},{child_aa},{mutation_category},{score},{prediction},{homoplasy},{occurrence}\n")                                    
        
        if not site_found:
            for snp_site in branch_snps_dict[site]:
                parent,child,snp,dimer = snp_site
                apobec = "False"
                if snp in ["C->T","G->A"] and dimer in ["GA","TC"]:
                    apobec = "True"
                fw.write(f"{site},NA,NA,{snp},{dimer},{apobec},NA,{parent},NA,NA,{child},NA,NA,intergenic,NA,NA,{homoplasy},{occurrence}\n")                                    
            
    fw.close()
            
        
def get_reconstruction_amino_acids(directory,alignment,treefile):
    node_states = get_node_states_all_sites(directory,f"{alignment}.state",alignment)
    print("Node states",len(node_states))
    branch_snps = f"{directory}/{treefile}.branch_snps.reconstruction.csv"

    reconstruct_amino_acid_mutations(f"{branch_snps}",
                                    node_states, 
                                     f"{directory}/{treefile}.amino_acid.reconstruction.csv")
    
    
import datetime
def year_fraction(date):
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    return date.year + float(date.toordinal() - start) / year_length

def get_root_to_tip_counts(aa_reconstruction,state_diffs,APO_out,root_to_tip_counts):

    site_info = {}
    with open(aa_reconstruction,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_info[row["site"]] = row

    # for i in site_info:
    #     print(i, site_info[i])

    seq_snps = collections.defaultdict(list)
    apobec_snps = collections.defaultdict(list)
    non_apo = collections.defaultdict(list)
    with open(state_diffs,"r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            root_variant = row["Node1"]

            for seq in reader.fieldnames:
                if seq != "site" and not seq.startswith("Node"):
                    if row[seq] != root_variant:
                        if "" not in [row[seq],root_variant]:
                            try:
                                snp_info = site_info[row["site"]]
                                if snp_info["apobec"] == "True":
                                    apobec_snps[seq].append(row["site"])
                                else:
                                    non_apo[seq].append(row["site"])
                                seq_snps[seq].append(row["site"])
                            except:
                                pass
    s = []
    date_apo = {}
    fw2 = open(APO_out,"w")
    fw2.write("name,date\n")
    with open(root_to_tip_counts,"w") as fw:
        fw.write("name,all_snps,apobec_snps,non_apobec_snps,date,decimal_year,precision\n")
        for i in seq_snps:
            datestring = i.split("|")[-1]
            
            precision = ""
            
            if len(datestring.split("-")) == 3:
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "day"
                print(datestring,odate,precision)
            elif len(datestring.split("-")) == 2:
                if datestring.endswith("-2"):
                    datestring+="-15"
                else:
                    datestring+="-16"
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "month"
                print(datestring,odate,precision)
            else:
                
                odate = float(datestring)+0.5
                precision="year"
                print(datestring,odate,precision)
            if i in apobec_snps:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},0,{datestring},{odate},{precision}\n")
                date_apo[i] = [apobec_snps[i],datestring]
            else:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},0,{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},0,0,{datestring},{odate},{precision}\n")

#             else:
#                 print(date)


def get_root_to_tip_counts_date_in(aa_reconstruction,state_diffs,APO_out,root_to_tip_counts,date_dict):

    site_info = {}
    with open(aa_reconstruction,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            site_info[row["site"]] = row

    # for i in site_info:
    #     print(i, site_info[i])

    seq_snps = collections.defaultdict(list)
    apobec_snps = collections.defaultdict(list)
    non_apo = collections.defaultdict(list)
    with open(state_diffs,"r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            root_variant = row["Node1"]

            for seq in reader.fieldnames:
                if seq != "site" and not seq.startswith("Node"):
                    if row[seq] != root_variant:
                        if "" not in [row[seq],root_variant]:
                            try:
                                snp_info = site_info[row["site"]]
                                if snp_info["apobec"] == "True":
                                    apobec_snps[seq].append(row["site"])
                                else:
                                    non_apo[seq].append(row["site"])
                                seq_snps[seq].append(row["site"])
                            except:
                                pass
    s = []
    date_apo = {}
    fw2 = open(APO_out,"w")
    fw2.write("name,date\n")
    with open(root_to_tip_counts,"w") as fw:
        fw.write("name,all_snps,apobec_snps,non_apobec_snps,date,decimal_year,precision\n")
        for i in seq_snps:
            datestring = date_dict[i]
            
            precision = ""
            
            if len(datestring.split("-")) == 3:
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "day"
                print(datestring,odate,precision)
            elif len(datestring.split("-")) == 2:
                if datestring.endswith("-2"):
                    datestring+="-15"
                else:
                    datestring+="-16"
                ddate = date.fromisoformat(datestring)
                odate = year_fraction(ddate)
                
                precision = "month"
                print(datestring,odate,precision)
            else:
                
                odate = float(datestring)+0.5
                precision="year"
                print(datestring,odate,precision)
            if i in apobec_snps:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},{len(apobec_snps[i])},0,{datestring},{odate},{precision}\n")
                date_apo[i] = [apobec_snps[i],datestring]
            else:
                if i in non_apo:
                    fw.write(f"{i},{len(seq_snps[i])},0,{len(non_apo[i])},{datestring},{odate},{precision}\n")
                else:
                    fw.write(f"{i},{len(seq_snps[i])},0,0,{datestring},{odate},{precision}\n")

def get_heptamers_internal(node_states,branch_snps):
    site_dict = collections.defaultdict(list)
    snp_count = 0
    apo_snp_count =0
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            
            if row["child"].startswith("Node") and row["parent"] != "Node1":
                snp_count +=1
                if row["snp"] in ["G->A","C->T"]:
                    apo_snp_count +=1
                    site_dict[row["site"]] = (row["parent"],row["snp"])
                    
    print(snp_count, apo_snp_count)
    heptamers = []
    for site in site_dict:
        
        parent,dimer = site_dict[site]
#         print(site, parent, dimer)
        site = int(site)
        indexes = []
        bases = []
        for i in range(site-3,site+4):
            indexes.append(i)
            nodes = node_states[f"{i}"]
            base = ""
            for node in nodes:
                if node[0] == parent:
                    base = node[1]
            bases.append(base)
        seq = "".join(bases)
        
        seq = Seq(seq)
        
        if dimer == "G->A":
            
            seq = seq.reverse_complement()
        
        heptamers.append(seq)
    return snp_count,heptamers
            
def get_heptamers_all(node_states,branch_snps):
    site_dict = collections.defaultdict(list)
    snp_count = 0
    with open(branch_snps,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            snp_count +=1
            if row["snp"] in ["G->A","C->T"]:
                site_dict[row["site"]] = (row["parent"],row["snp"])
    heptamers = []
    for site in site_dict:
        
        parent,dimer = site_dict[site]
#         print(site, parent, dimer)
        site = int(site)
        indexes = []
        bases = []
        for i in range(site-3,site+4):
            indexes.append(i)
            nodes = node_states[f"{i}"]
            base = ""
            for node in nodes:
                if node[0] == parent:
                    base = node[1]
            bases.append(base)
        seq = "".join(bases)
        
        seq = Seq(seq)
        
        if dimer == "G->A":
            
            seq = seq.reverse_complement()
        
        heptamers.append(seq)
    return snp_count,heptamers

def make_heptamer_figure(heptamers,hept_out,heptfigure):
    hept_dict = {}
    for i in range(7):
        hept_dict[i] = collections.Counter()

    rel_dict = {0:-3,1:-2,2:-1,3:0,4:1,5:2,6:3}
    mutations = 0

    with open(hept_out,"w") as fw:
        fw.write("position,rel_position,base\n")
        for heptamer in heptamers:
            hept_list = list(heptamer)

            for i in range(len(hept_list)):
                hept_dict[i][hept_list[i]] +=1 
                fw.write(f"{i+1},{rel_dict[i]},{hept_list[i]}\n")
                
#             for i in hept_dict:
#                 print(i, hept_dict[i])
                        
    df = pd.read_csv(hept_out)
    df
    fig,ax= plt.subplots(figsize=(4,4),facecolor='w',frameon=False)
    colors = {"A":"#4682B4",
              "C":"indianred",
              "G":"#BC9D60",
              "T":"#87CEEB"}
    
    colors = ["#BC9D60","#87CEEB","#4682B4","indianred"]
    customPalette = sns.set_palette(sns.color_palette(colors))

#     legend_elements = [Line2D([0], [0], markerfacecolor='#4682B4',color='white',markersize=15,alpha=0.6, marker='o',label='A'),
#                        Line2D([0], [0], markerfacecolor='indianred',color='white',markersize=15,alpha=0.6, marker='o',label='C'),
#                        Line2D([0], [0], markerfacecolor='#BC9D60',color='white',markersize=15, alpha=0.6,marker='o',label='G'),
#                         Line2D([0], [0], markerfacecolor='#87CEEB',color='white',markersize=15,alpha=0.6, marker='o',label='T')]

    g = sns.histplot(data=df, x="rel_position",hue="base",multiple="fill",hue_order = ['A','C','G','T'],palette=customPalette,alpha=0.6,bins=7,linewidth=1)

    sns.move_legend(ax, "center", bbox_to_anchor=(1, 1), frameon=False,title="")
    plt.xlim(-3,3)

    plt.xlabel("Position\n");
    plt.ylabel("Proportion")
    [ax.spines[loc].set_visible(False) for loc in ['top','right']]

    plt.savefig(f"figures/{heptfigure}.svg");
    plt.savefig(f"figures/{heptfigure}.png",bbox_inches='tight', 
                   transparent=True);
    plt.show();
    
def get_reconstruction_heptamers(directory, alignment,treefile,hept_out,heptamer_figure,width=35):
    node_states = get_node_states_all_sites(directory,f"{alignment}.state",alignment)
    print("Node states",len(node_states))
    branch_snps = f"{directory}/{treefile}.branch_snps.reconstruction.csv"
    snp_count,heptamers = get_heptamers_all(node_states,branch_snps)

#     snp_count,heptamers = get_heptamers_internal(node_states,branch_snps)
    print("SNP count:",snp_count)
    print("Heptamer count:",len(heptamers))
    hept_out = f"{directory}/{hept_out}"
    hept_figure = f"{heptamer_figure}"

    make_heptamer_figure(heptamers,hept_out,hept_figure)
    
    

def get_reconstruction_internal_heptamers(directory, alignment,treefile,hept_out,heptamer_figure,width=35):
    node_states = get_node_states_all_sites(directory,f"{alignment}.state",alignment)
    print("Node states",len(node_states))
    branch_snps = f"{directory}/{treefile}.branch_snps.reconstruction.csv"
#     snp_count,heptamers = get_heptamers_all(node_states,branch_snps)

    snp_count,heptamers = get_heptamers_internal(node_states,branch_snps)
    print("SNP count:",snp_count)
    print("Heptamer count:",len(heptamers))
    hept_out = f"{directory}/{hept_out}"
    hept_figure = f"{heptamer_figure}"
    make_heptamer_figure(heptamers,hept_out,hept_figure)



def make_non_syn_tree_fig(outfile,branch_snps,treefile,width,height):
    
    branch_snps_dict = read_in_branch_snps(branch_snps)
    print(len(branch_snps_dict),"branches")
#     acc_dict = get_acc_to_metadata_map(metadata)
#     print(len(acc_dict))        
    
    my_tree=bt.loadNewick(treefile,absoluteTime=False)

    fig,ax = plt.subplots(figsize=(width,height),facecolor='w')

    x_attr=lambda k: k.height ## x coordinate of branches will be absoluteTime attribute
    su_func=lambda k: 50-30*k.height/my_tree.treeHeight ## size of tips
    s_func=lambda k: 50-20*k.height/my_tree.treeHeight ## size of tips
    c_func=lambda k: "dimgrey"

    increment = my_tree.treeHeight/150
    my_tree.plotTree(ax,x_attr=x_attr) ## plot branches
    my_tree.plotPoints(ax,size=s_func,colour=c_func,x_attr=x_attr) ## plot circles at tips
    mpl.rcParams['font.family'] = 'sans-serif'

    for k in my_tree.Objects:
        current_node = k
        if k.branchType == 'leaf':
            current_node.traits["label"]=k.name

        node_name = current_node.traits["label"]
        try:
            parent_name = current_node.parent.traits["label"]
        except:
            continue
        branch_name= f"{parent_name}_{node_name}"

        if branch_name in branch_snps_dict:
            snps = []
#                 print(branch_name, len(branch_snps_dict[branch_name]))
            snp_placement = current_node.parent.height + increment
            for s in branch_snps_dict[branch_name]:
                site,snp,dimer = s
                if snp == "G->A":
                    if dimer in ["GA"]:
                        snps.append((1,"orange"))
                    else:
                        snps.append((3,"lightseagreen"))
                elif snp == "C->T":
                    if dimer in ["TC"]:
                        snps.append((2,"blue"))
                    else:
                        snps.append((3,"lightseagreen"))
                else:
                    snps.append((4,"dimgrey"))

            for snp in sorted(snps, key = lambda x : x[0]):
                plt.scatter([snp_placement],[k.y+0.5],color=snp[1],s=30)
                snp_placement += increment

    [ax.spines[loc].set_visible(False) for loc in ['top','right','left','bottom']]
    ax.tick_params(axis='y',size=0)
    ax.tick_params(axis='x',size=0)

    ax.set_yticklabels([])
    ax.set_xticklabels([])

    plt.savefig(f"figures/{outfile}.svg");
    plt.savefig(f"figures/{outfile}.png",bbox_inches='tight', 
                   transparent=True);
    plt.show()


def forward_amino_acid_mutations(type_site_dict,seq,type_counter):
    genes = get_gene_boundaries()
    grantham_scores = get_grantham_scores()
    
    
    for i in range(len(seq)):
        base_1 = seq[i-1]
        base_2 = seq[i]
        
        dimer = f"{base_1}{base_2}"
        edited = "false"
        if dimer in ["TC","GA"]:
            snp = ""
            site = ""
            parent = ""
            child = ""
            if dimer == "TC":
                snp = "C->T"
                site = i
                parent = "C"
                child = "T"
            else:
                snp = "G->A"
                site = i-1
                parent = "G"
                child = "A"
            
            site_found = False
            
            for gene in genes:
                start,end=gene
                if site_found == True:
                    break
                    
                if site in range(start,end):
                    
                    name,length,direction = genes[gene]
                    if direction == "forward":
                        aa_position = get_aa_position(site,gene[0])
                        codon_indexes = get_codon_indexes(aa_position,site)
                    else:
                        aa_position = reverse_aa_position(start,end,site)
                        codon_indexes = get_codon_indexes_rev_strand(aa_position,site)

                    parent_codon = []
                    child_codon = []
                    for i in codon_indexes:
                        if i == site:
                            child_codon.append(child)
                        else:
                            child_codon.append(seq[i])
                        parent_codon.append(seq[i])


                    parent_codon = "".join(parent_codon)
                    child_codon = "".join(child_codon)

                    parent_codon = Seq(parent_codon)
                    child_codon = Seq(child_codon)
                    
                    if direction == "reverse":
                        parent_codon= parent_codon.reverse_complement()
                        child_codon= child_codon.reverse_complement()

                    parent_aa = parent_codon.translate()
                    child_aa = child_codon.translate()
                    mutation_category,score,prediction = categorise_amino_acid_mutation(parent_aa,child_aa,grantham_scores)
                    
                    type_counter[mutation_category] +=1
                    type_site_dict[mutation_category].append(site)
                    site_found = True

            if not site_found:
                mutation_category="intergenic"
                type_counter[mutation_category] +=1
                type_site_dict[mutation_category].append(site)
            

def run_full_analysis(directory, alignment, tree):
    generate_reconstruction_files(directory, 
                                  alignment,
                                  tree)

    load_info(directory,alignment,tree,f"{tree}.fig")

    make_apobec_context_internal_plot(f"{directory}/{tree}.branch_snps.reconstruction.csv",
                                            f"{directory}/{tree}.branch_snps.internal.counts.csv",
                                            f"{tree}.internal.counts")


    make_apobec_context_mutation_count_plot(f"{directory}/{tree}.branch_snps.reconstruction.csv",
                                            f"{directory}/{tree}.branch_snps.counts.csv",
                                            f"{tree}.all.counts")

    make_apobec_context_external_plot(f"{directory}/{tree}.branch_snps.reconstruction.csv",
                                            f"{directory}/{tree}.branch_snps.external.counts.csv",
                                            f"{tree}.external.counts")

    get_reconstruction_amino_acids(directory,alignment,tree)

    get_root_to_tip_counts(f"{directory}/{tree}.amino_acid.reconstruction.csv",
                           f"{directory}/{tree}.state_differences.csv",
                           f"{directory}/APOBEC_reconstructed_SNPs.csv",
                           f"{directory}/root_to_tip.data.csv")



my_ref_file = "/Users/s1680070/repositories/alignHPXV/squirrel/data/NC_063383.aln.fasta"

directory = "updated_2023-02-03/"
tree = "Clade_IIb_2023-02-03_with_gisaid.aln.fasta.pruned.tree"
alignment="Clade_IIb_2023-02-03_with_gisaid.aln.fasta"

run_full_analysis(directory, alignment, tree)
