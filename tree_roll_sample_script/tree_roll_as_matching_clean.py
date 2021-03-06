##########################################
#Author: Ajay
#Description: Probabilistic matching of single cell gene expression vectors to appropriate Anatomical Structures guided by EMAP-GXD DAG Tree
#Input Arguments: theiler stage e.g. "TS13"

##########################################

import networkx as nx
import sys
import numpy as np
import os
from copy import deepcopy
from collections import Counter
import subprocess
import re

#Function to store gene expression in form of dictionary.
def processGeneTS13(fileName, ts):
    detected = {'Moderate': 'Moderate', 'Trace': 'Trace',
                'Weak': 'Weak', 'Very strong': 'Very strong', 'Strong': 'Strong',
                'Present': 'Present'}

    notDetected = {'Absent': 'Absent'}

    tissueExpression = {}

    with open(fileName, 'r') as inp:

        inp.readline()  #skip first line

        for line in map(str.strip, inp.readlines()):

            values = line.split('\t')
            if values[1] == '':  #choose non-mutant data

                if values[3] == ts[2:]:

                    if values[5] in detected or values[5] in notDetected:

                        if values[4] not in tissueExpression:
                            tissueExpression[values[4]] = {}  #tissueExpression['tissue']={}
                        if values[0] not in tissueExpression[values[4]]:
                            tissueExpression[values[4]][values[0]]=[]
                        if values[5] in detected:
                            tissueExpression[values[4]][values[0]].append(1)  #tissueExpression['tissue']['gene']={}
                        if values[5] in notDetected:
                            tissueExpression[values[4]][values[0]].append(0)

    for AS in tissueExpression:
        for gene in tissueExpression[AS]:
            if len(tissueExpression[AS][gene])>1:


                data = Counter(tissueExpression[AS][gene]) #Function to cound frequency of different elements in a list

                count_data=data.most_common(2)  # Returns list of tupules with frequency of 1 and 0 [(0,#),(1,#)]
                if len(count_data)==2:
                    if count_data[0][1]==count_data[1][1]:
                        tissueExpression[AS][gene]=[1]
                    else:
                        tissueExpression[AS][gene]=[data.most_common(1)[0][0]]
                else:
                    tissueExpression[AS][gene]=[tissueExpression[AS][gene][0]]



    return tissueExpression

#Function to build DAG tree and annotate nodes with gene expression info. 
def buildGraph(fileName, tissueExpression, theiler_stage):
    directed_graph = nx.DiGraph()

    startTerm = '[Term]'

    lines = []
    with open(fileName, 'r') as inp:
        lines = map(str.strip, inp.readlines())

    linePos = 0
    countTS = 0
    TSInfo = {}
    edges = []
    while (linePos < len(lines)):

        line = lines[linePos]

        if startTerm in line:
            linePos += 1
            emapID = lines[linePos].split(' ')[1].strip()
            emapID = emapID.replace(':', '_')
            linePos += 1
            name = lines[linePos].split(':')

            if theiler_stage in name[1]:

                tissue = ' '.join(name[1].strip().split(' ')[1:])
                #tissue = ' '.join(name[1:]).strip()
                #print '%s - %s' % (emapID, tissue)
                countTS += 1

                if emapID not in TSInfo:
                    tissue = tissue.strip()
                    emapID = emapID.strip()
                    TSInfo[emapID] = tissue
                    attr_dict = {}
                    attr_dict['name'] = tissue
                    expression_num = 0
                    if tissue in tissueExpression:
                        expression_num = len(tissueExpression[tissue].keys())
                    attr_dict['expression_number'] = expression_num

                    directed_graph.add_node(emapID, attr_dict)

                contRead = True
                while ( contRead ):
                    linePos += 1
                    if 'relationship' in lines[linePos] and 'part_of' in lines[linePos]:
                        partOf = lines[linePos].split('part_of')[1].split('!')[0].strip()
                        #relationship[partOf] = emapID
                        #print '%s\t%s' % (partOf, emapID)
                        partOf = partOf.replace(':', '_')
                        edges.append((partOf, emapID))
                    elif 'is_a' in lines[linePos]:
                        pass
                    else:
                        contRead = False
        linePos += 1

        directed_graph.add_edges_from(edges)

    return directed_graph

#Function to prune DAG tree based on minimum number of markers
def prune_graph(directed_graph, root_node, min_number):
    check = []
    temp = nx.bfs_edges(directed_graph, root_node)
    for test in temp:

        check.append(test)

    check.reverse()

    for test in check:

        parent = test[0]
        child = test[1]




        if directed_graph.node[child]['expression_number'] < min_number:
            for successor in directed_graph.successors(child):
                directed_graph.remove_edge(child, successor)
                directed_graph.add_edge(parent, successor)
            directed_graph.remove_edge(parent, child)
            directed_graph.remove_node(child)

    return directed_graph

#Function to merge two dictionaries together.
def merge_dols(dol1, dol2):
  keys = set(dol1).union(dol2)
  no = []
  return dict((k, dol1.get(k, no) + dol2.get(k, no)) for k in keys)

#Function to impute gene expression in each node as union of the nodes below it (in the trajectory).
def collapse_exp_node_new2(directed_graph, root_node, expression):
    check = []
    directed_graph.node[root_node]['count_cells']=0
    directed_graph.node[root_node]['gene_collapsed']={}
    temp = nx.bfs_edges(directed_graph, root_node)
    for test in temp:

        check.append(test)

    check.reverse()

    for test in check:
        parent = test[0]

        child = test[1]



        child_as=directed_graph.node[child]['name']

        genes_child=deepcopy(expression[child_as])
        directed_graph.node[child]['count_cells']=0
        if child_as=='head' or child_as=='trunk':
            #directed_graph.node[child]['count_cells']=0
            directed_graph.node[child]['gene_collapsed']=deepcopy(genes_child)
            continue


        if len(directed_graph.successors(child))>0:              #If child has descendants
            for successor in directed_graph.successors(child):
                as_name=directed_graph.node[successor]['name']

                genes_successor=directed_graph.node[successor]['gene_collapsed']
                genes_child=merge_dols(genes_child,genes_successor)

        if 'gene_collapsed' not in directed_graph.node[child].keys():    #if no instance of collapsed gene annotation

            directed_graph.node[child]['gene_collapsed']=deepcopy(genes_child)
            directed_graph.node[child]['count_cells']=0
        else:
            directed_graph.node[child]['gene_collapsed']=merge_dols(directed_graph.node[child]['gene_collapsed'],genes_child)

    check.reverse()

    return directed_graph

#Print out file containing pairwise lineage relationships
def lineage_file_print(directed_graph, root, lineageFile):
    lineageOut = open(lineageFile, 'w')

    edgeList = nx.bfs_edges(directed_graph, root)

    for oParent, oChild in edgeList:  #directed_graph.edges_iter():
        #print directed_graph.node[oParent]['name'],directed_graph.node[oChild]['name']
        try:
            parentName = directed_graph.node[oParent]['name'] + '_' + str(len(directed_graph.node[oParent]['gene_collapsed'])) + '_(' + str(directed_graph.node[oParent]['count_cells'])+')'
            childName = directed_graph.node[oChild]['name'] + '_' + str(len(directed_graph.node[oChild]['gene_collapsed'])) + '_(' + str(directed_graph.node[oChild]['count_cells'])+')'
            lineageOut.write('%s\t%s\n' % (parentName, childName))
        except KeyError:
            #print directed_graph.node[oParent]['name']
            parentName = directed_graph.node[oParent]['name']+ '_' + str(len(directed_graph.node[oChild]['expression_number']))+'_(' + str(directed_graph.node[oChild]['count_cells'])+')'

            if directed_graph.node[oChild]['name'] == 'head' or directed_graph.node[oChild]['name']== 'trunk':
                childName = directed_graph.node[oChild]['name'] +'_' + str(len(directed_graph.node[oChild]['expression_number'])) +'_(' + str(directed_graph.node[oChild]['count_cells'])+')'
            else:
                childName = directed_graph.node[oChild]['name'] + '_' + str(len(directed_graph.node[oChild]['gene_collapsed'])) + '_(' + str(directed_graph.node[oChild]['count_cells'])+')'
            lineageOut.write('%s\t%s\n' % (parentName, childName))
            continue
    lineageOut.close()

#Output expression table file from modified graph object
def output_collapsed_expression_tsv (directed_graph,root,out,expression):
    out=open(out,'w')
    genes_96='Aldh1a2,Cdx4,Dlx2,Fgf15,Fut4,Gtf2ird1,Hoxb7,Kit,Nog,Prtg,Smad2,Tcf15,Bmp2,Cgnl1,Dmbx1,Fgf4,Gata4,Hey1,Hoxb9,Lhx5,Olfm1,Rarb,Snai1,Tfap2a,Bmp4,Crabp1,Egr2,Fgf8,Gbx1,Hey2,Hoxd1,Mafb,Otx2,Robo3,Snai2,Tll1,Bmp7,Crabp2,Epha2,Flt4,Gbx2,Hhex,Irx3,Meox1,Pax6,Scube2,Sox1,Tmed2,Bmpr1b,Cyp26c1,Epha7,Foxa2,Gpc4,Hoxa1,Itga4,Mtf2,Pdgfra,Sdc1,Sox2,Twist1,Cdh1,Dach1,Ephb2,Foxd3,Grhl1,Hoxa3,Itga6,Myl2,Pecam1,Shh,Sox9,Wnt1,Cdh11,Dll1,Ephb3,Foxf1,Grhl2,Hoxb1,Kdm5b,Ncam1,Prrx1,Sik1,T,Wnt3a,Cdon,Dll3,Evx1,Fst,Grhl3,Hoxb2,Kdr,Nodal,Prrx2,Six1,Tbx6,Zmiz2'
    gene_list=genes_96.rsplit(',') #Make gene list

    out.write(''+'\t'+"\t".join(gene_list)+'\n')
    edgeList = nx.bfs_edges(directed_graph, root)

    list_done=[]
    for oParent, oChild in edgeList:  #directed_graph.edges_iter():
        parentName = directed_graph.node[oChild]['name']
        if parentName in list_done:
            continue

        list_done.append(parentName)
        print_gene_exp=[]

        print_gene_exp.append(oChild) # prints AS_id
        if parentName=='head'or parentName=='trunk':
            for gene in gene_list:
                try:
                    if 1 in expression[parentName][gene]:
                        print_gene_exp.append('1')
                    else:
                        print_gene_exp.append('0')
                except KeyError:
                    print_gene_exp.append('NA')
            print_line="\t".join(print_gene_exp)
            out.write(print_line+'\n')
            continue
        for gene in gene_list:
            try:
                if 1 in directed_graph.node[oChild]['gene_collapsed'][gene]:
                    print_gene_exp.append('1')
                else:
                    print_gene_exp.append('0')
            except KeyError:
                print_gene_exp.append('NA')
        print_line="\t".join(print_gene_exp)
        out.write(print_line+'\n')   # indent ends here.

#print in order
def list_to_qpcr96(gene_list,gene_dict):

    print_gene_exp=[]
    for gene in gene_list:
        try:
            if 1 in gene_dict[gene]:
                print_gene_exp.append('1')
            else:
                print_gene_exp.append('0')
        except KeyError:
            print_gene_exp.append('NA')
    print_line="\t".join(print_gene_exp)
    return print_line

#python wrapper function to R script
def euc_calc_R (directed_graph,root,expression,single_cell_gene_list):

    temp_file='output/temp_table.csv'
    temp_out=open(temp_file,'w')
    genes_96='Aldh1a2,Cdx4,Dlx2,Fgf15,Fut4,Gtf2ird1,Hoxb7,Kit,Nog,Prtg,Smad2,Tcf15,Bmp2,Cgnl1,Dmbx1,Fgf4,Gata4,Hey1,Hoxb9,Lhx5,Olfm1,Rarb,Snai1,Tfap2a,Bmp4,Crabp1,Egr2,Fgf8,Gbx1,Hey2,Hoxd1,Mafb,Otx2,Robo3,Snai2,Tll1,Bmp7,Crabp2,Epha2,Flt4,Gbx2,Hhex,Irx3,Meox1,Pax6,Scube2,Sox1,Tmed2,Bmpr1b,Cyp26c1,Epha7,Foxa2,Gpc4,Hoxa1,Itga4,Mtf2,Pdgfra,Sdc1,Sox2,Twist1,Cdh1,Dach1,Ephb2,Foxd3,Grhl1,Hoxa3,Itga6,Myl2,Pecam1,Shh,Sox9,Wnt1,Cdh11,Dll1,Ephb3,Foxf1,Grhl2,Hoxb1,Kdm5b,Ncam1,Prrx1,Sik1,T,Wnt3a,Cdon,Dll3,Evx1,Fst,Grhl3,Hoxb2,Kdr,Nodal,Prrx2,Six1,Tbx6,Zmiz2'
    gene_list=genes_96.rsplit(',') #Make gene list
    temp_out.write(''+'\t'+"\t".join(gene_list)+'\n')
    temp_out.write("\t".join(single_cell_gene_list)+'\n')
    if len(directed_graph.successors(root))==1:

        gate="on"
        as_node=directed_graph.successors(root)
        as_node=as_node[0]
        directed_graph.node[as_node]['count_cells']+=1

        return directed_graph,gate,as_node


    if len(directed_graph.successors(root))>1: #if node has descendants
        descend = directed_graph.successors(root)
        for x in descend:
            AS_name=directed_graph.node[x]['name']

            try:
                as_expression=directed_graph.node[x]['gene_collapsed']
            except KeyError:
                as_expression= expression[AS_name]
            print_line=list_to_qpcr96(gene_list,as_expression)
            temp_out.write(x+'\t'+print_line+'\n')

        temp_out.close()
        back_data='output/back.rds'
        cmd= 'Rscript tree_matching.R {0}/{1} {2}/{3}'.format(curr_dir, back_data, curr_dir, temp_file)

        result = subprocess.check_output(cmd, shell=True)
        as_node=re.findall('"([^"]*)"',result)
        as_node=as_node[0]
        gate="on"

        directed_graph.node[as_node]['count_cells']+=1
        return directed_graph,gate,as_node


    else:

        gate='off'

        return directed_graph,gate,root

#Tree roll function
def tree_roll_master(cell_data,directed_graph,root_embryo,expression,tree_roll_exp_annot_out):
    count_cell=0
    with open(cell_data, 'r') as inp, open(tree_roll_exp_annot_out, 'w') as out_exp:

        first_line=inp.readline()  #skip first line
        out_exp.write(first_line)
        for line in map(str.strip, inp.readlines()):
            count_cell+=1
            gate='on'
            values = line.split(',')
            single_cell_gene_list=values[0:97]
            as_node=root_embryo
            while gate == "on":
                directed_graph,gate,as_node=euc_calc_R(directed_graph,as_node,expression,single_cell_gene_list)
                values.append(directed_graph.node[as_node]['name'])
            #print count_cell
            out_exp.write(','.join(values)+'\n')



    return directed_graph

# Define Legend for the output graph
def assignColorKeys(directed_graph, max_limit):

    ranges = {}
    ranges[0] = []
    ranges[1] = []
    ranges[2] = []
    ranges[3] = []
    ranges[4] = []
    ranges[5] = []
    ranges[6] = []
    ranges[7] = []
    ranges[8] = []
    ranges[9] = []

    cellNumColor = {}

    #matchColors = ['#0000CD', '#0066FF', '#3385FF', '#87CEFF', '#EBF5FF', '#FFEBE6', '#FF704D', '#FF4719', '#E62E00', '#B22400' ]
    matchColors = [  "#0000FF", "#1C00E2", "#3800C6", "#5500AA", "#71008D", "#8D0071", "#AA0055", "#C60038", "#E2001C", "#FF0000" ]
    #matchColors = ['#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF','#FFFFFF']


    #numCells=np.array( [directed_graph.node[node]['gene_collapsed'] for node in directed_graph.node], dtype=float )

    numCells=np.array([],dtype=float)
    for node in directed_graph.node.keys():
        try:
            numCells=np.append(numCells,directed_graph.node[node]['count_cells'])

        except KeyError:
            continue


    divideby = max_limit - numCells.min()
    countGT = 0

    for value in numCells:

        if value > max_limit:
            countGT += 1
            continue

        scale = float( (value - numCells.min()) / divideby )

        if scale >= 0.1 and scale < 0.2:
            color = matchColors[1]
            ranges[1].append(value)
        elif scale >= 0.2 and scale < 0.3:
            color = matchColors[2]
            ranges[2].append(value)
        elif scale >= 0.3 and scale < 0.4:
            color = matchColors[3]
            ranges[3].append(value)
        elif scale >= 0.4 and scale < 0.5:
            color = matchColors[4]
            ranges[4].append(value)
        elif scale >= 0.5 and scale < 0.6:
            color = matchColors[5]
            ranges[5].append(value)
        elif scale >= 0.6 and scale < 0.7:
            color = matchColors[6]
            ranges[6].append(value)
        elif scale >= 0.7 and scale < 0.8:
            color = matchColors[7]
            ranges[7].append(value)
        elif scale >= 0.8 and scale < 0.9:
            color = matchColors[8]
            ranges[8].append(value)
        elif scale >= 0.9 and scale <= 1.0:
            color = matchColors[9]
            ranges[9].append(value)
        else:
            color = matchColors[0]
            ranges[0].append(value)

        if value not in cellNumColor:
            cellNumColor[value] = color



    legend2print = []
    legend2print.append( 'Legend [shape=none, margin=0, label=<' )
    legend2print.append( '<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">' )
    legend2print.append( '<TR>' )
    legend2print.append( ' <TD COLSPAN="2"><B>Legend</B></TD>' )
    legend2print.append( '</TR>' )

    for x in xrange(len(ranges)):

        if len(ranges[x]) == 0: continue
        addLegend = '<TR><TD>%s - %s Cells</TD> <TD BGCOLOR="%s"> %s Anatomical Struct(s) </TD> </TR>' % (min(ranges[x]) , max(ranges[x]), matchColors[x], len(ranges[x]))
        legend2print.append( addLegend )

    addLegend = '<TR><TD> More than %s Cells</TD> <TD BGCOLOR="%s"> %s Anatomical Struct(s) </TD> </TR>' % (max_limit , 'green', countGT)
    legend2print.append( addLegend )

    legend2print.append( '</TABLE>' )
    legend2print.append( '>];' )

    return cellNumColor, legend2print

# Output graph file
def output_graph(cell_graph, root, outFile, cellColor, legend2print):


    out = open(outFile, 'w')

    white = '#ffffff'
    fontSize = '10'

    edgeList = nx.bfs_edges(cell_graph, root)

    out.write( 'digraph test123 {\n rankdir=LR;\n subgraph cluster01 { \n  ' )
    nodes = {}
    for oParent, oChild in edgeList: #directed_graph.edges_iter():

        parentName = cell_graph.node[oParent]['name']
        childName = cell_graph.node[oChild]['name']


        parent = oParent
        if parent not in nodes:

            nodes[parent] = {}

            tempLabel = cell_graph.node[oParent]['name']
            if len(tempLabel) > 10 and ' ' in tempLabel:
                temp = tempLabel.split(' ')
                if len(temp) > 3:
                    newLabel = '%s %s\\n%s' % (temp[0], temp[1], ' '.join(temp[2:]))
                else:
                    newLabel = '%s\\n%s' % (temp[0], ' '.join(temp[1:]))

            else:
                newLabel = tempLabel

            nodes[parent]['label'] = newLabel
            nodes[parent]['gene_collapsed'] = len(cell_graph.node[oParent]['gene_collapsed'])
            nodes[parent]['count_cells']= cell_graph.node[oParent]['count_cells']


        child = oChild
        if child not in nodes:

            nodes[child] = {}

            tempLabel = cell_graph.node[oChild]['name']
            if len(tempLabel) > 10 and ' ' in tempLabel:
                temp = tempLabel.split(' ')
                if len(temp) > 3:
                    newLabel = '%s %s\\n%s' % (temp[0], temp[1], ' '.join(temp[2:]))
                else:
                    newLabel = '%s\\n%s' % (temp[0], ' '.join(temp[1:]))
            else:
                newLabel = tempLabel

            nodes[child]['label'] = newLabel
            nodes[child]['gene_collapsed'] = len(cell_graph.node[oChild]['gene_collapsed'])
            nodes[child]['count_cells'] = cell_graph.node[oChild]['count_cells']

        out.write( '%s -> %s;\n' % ( parent, child ) )

    for node in nodes:

        number = ''
        mark = ''
        colorCode = white
        if nodes[node]['count_cells'] > 0:
            number = '_%s_(%s)' % (nodes[node]['gene_collapsed'],nodes[node]['count_cells'])

            colorCode = 'green'
            if nodes[node]['count_cells'] in cellColor:
                colorCode = cellColor[ nodes[node]['count_cells'] ]
        else:
            number = '_%s_(%s)' % (nodes[node]['gene_collapsed'],nodes[node]['count_cells'])

        mark = ',style="filled",fillcolor="%s"' % (colorCode) #(geneNumColor[ nodes[node]['expression'] ])
        out.write( '%s[label="%s%s"%s];\n' % (node, nodes[node]['label'], number, mark) )

    out.write(  '}\n' )

    out.write('subgraph cluster02 {\n')
    for line in legend2print:
        out.write( '%s\n' % line )
    out.write('}\n')


    out.write(  '}\n' )
    out.close()

if __name__ == '__main__':

    ######Enter input files#######
    gxd_file ='reference_data/biomark_96_mart_export_03232015.tsv' #Reference gene expression table from GXD
    dag_file ='reference_data/EMAP_combined.obo' #DAG Tree flat file
    exp_col_file = 'output/TS13_exp_col_prune_1_as_id.tsv' #Output reference gene expression file after DAG tree manipulation
    cell_data='experiment_data/all_plates_c5to12_raw_ct_bool_sample.csv' #Single cell boolean gene expression data
    tree_roll_exp_annot_out='output/treeroll_assignments.tsv' #Output of anatomical structure assignment (in each level of DAG tree) for each single cell vector
    graph_out_file = 'output/graph_output_file.txt' #file containing instructions to construct annotated DAG tree
    graph_out_fig='output/graph_output_figure.png'
    curr_dir=os.getcwd()
    ##############################

    ts = sys.argv[1] # user input Theiler Stage
    tissueExpression = processGeneTS13(gxd_file, ts) #store gene expression data

    directed_graph = buildGraph(dag_file, tissueExpression, ts) # Construct DAG Tree

    root_embryo = None
    root_extraembryonic = None
    for nd in directed_graph.nodes():

        if directed_graph.node[nd]['name'] == 'embryo': 
            root_embryo = nd  #set root node as embryo

        if directed_graph.node[nd]['name'] == 'extraembryonic component':
            root_extraembryonic = nd #set another root node as extra embryonic

    new_directed_graph = prune_graph(directed_graph, root_embryo, 1) # Remove Anatomical structures with no gene expression data

    new_directed_graph2=collapse_exp_node_new2(new_directed_graph, root_embryo, tissueExpression) #Impute data from daughter nodes -> parent nodes

    output_collapsed_expression_tsv(new_directed_graph2,root_embryo,exp_col_file,tissueExpression) #Ouput reference gene expression table post DAG tree manipulation

    os.system('Rscript create_Robject.R '+curr_dir+'/'+cell_data+' '+curr_dir+'/'+exp_col_file) #Create and store background model from gene expression data

    cell_graph=tree_roll_master(cell_data,new_directed_graph2,root_embryo,tissueExpression,tree_roll_exp_annot_out) #perform matching with rolling the single cell vectors down the DAG tree.

    cellNumColor, legend2print = assignColorKeys(cell_graph, 30.0) #legends for graph

    output_graph(cell_graph, root_embryo, graph_out_file, cellNumColor, legend2print) #output graph file

    os.system('cd '+curr_dir+';dot -Tpng '+graph_out_file+' > '+graph_out_fig) #generate image of graph tree with annotation




