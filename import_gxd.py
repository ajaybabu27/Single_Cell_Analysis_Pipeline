#Author: Ajay
#Description: Obtains wild type gene expression data from mouse GXD database.
#Input Arguments: file path for full GXD data file to be written.

from intermine.webservice import Service
import sys
service = Service("http://www.mousemine.org/mousemine/service")
query = service.new_query("GXDExpression")
query.add_view(
    "age", "stage", "strength", "detected", "feature.symbol", "structure.name",
    "genotype.symbol"
)

out=open(sys.argv[1],'w')
out.write('Gene Symbol	Mutant Allele(s)	Age	Theiler Stage	Anatomical Structure	Level\n')

for row in query.rows():

    if row["genotype.symbol"] == None:
        mutant=''
    stage=row["stage"][2:]
    if stage[0]=='0':
        stage=stage[1]
    out.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(row["feature.symbol"],mutant,row["age"], stage,row["structure.name"], row["strength"]))

out.close()


