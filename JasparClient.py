#! /usr/bin/python
#
#   Python script that communicates with the JasparDB web service at CBU
#   WSDL: http://api.bioinfo.no/wsdl/JasparDB.wsdl
#   By Jan Christian Bryne (chrb@ii.uib.no).
#

import SOAPpy, sys
from optparse import OptionParser
from pprint import pprint

def main():
    """Command line interface to the JasparDB database Web service """

    # command line parameter parsing configuration
    usage = "usage: JasparClient.py -d database -f format   -i/-n/-s/-a/"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--identifier", action="store", dest="identifier", help="retrive a profile by identifier")
    parser.add_option("-n", "--name", action="store", dest="name", help="retrive a profile by name")
    parser.add_option("-s", "--search", action="store", dest="search", help="search for profile matching given tags. Ex: tag=value1,value2")
    parser.add_option("-a", "--all_profiles", action="store_true", default=False, help="retrive all profiles")
    parser.add_option("-d", "--database", action="store", dest="database", help="retrive profile from given database. Either CORE, CNE, FAM or PHYLOFACTS")
    parser.add_option("-f", "--format", action="store", dest="format", help="retrive profile in given format. Either PFM, PWM or ICM. PFM is default")
    parser.add_option("-v", dest='verbose', action="store_true", default=False, help="print details of Web service communication")
    (options, args) = parser.parse_args()

   # web service client configuration    
    endpoint = 'http://api.bioinfo.no/services/JasparDB'
    namespace = 'http://bccs.uib.no/Jaspar'
    jaspar = SOAPpy.SOAPProxy(endpoint)
    jaspar.namespace = namespace
    jaspar.noroot = 1
    if options.verbose:
        jaspar.config.dumpSOAPOut = 1
        jaspar.config.dumpSOAPIn = 1

    try:
        format = SOAPpy.Types.untypedType(options.format)
        format._name = 'ns1:Format'
        database = SOAPpy.Types.untypedType(options.database)
        database._name = 'ns1:Database'
        
        if options.name:
            name = SOAPpy.Types.untypedType(options.name)
            name._name = 'ns1:Name'
            matrix = jaspar.getMatrixByName(name, format, database)
	    _print_matrix(matrix)
        
        if options.identifier:
            identifier = SOAPpy.Types.untypedType(options.identifier)
            identifier._name = 'ns1:Id'
            matrix = jaspar.getMatrixById(identifier, format, database)
            _print_matrix(matrix)
                        
        if options.all_profiles:
            matrices = jaspar.getAllMatrices(format, database)
            for matrix in matrices:
                _print_matrix(matrix)
        
        if options.search:
            query_string = options.search
            query = "";
            query += """<ns1:Name>%s</ns1:Name>""" % query_string.split('=')[0]
            values_line = query_string.split('=')[1]
            for value in values_line.split(','):
		      query += """<ns1:Value>%s</ns1:Value>""" % value
            query = SOAPpy.Types.untypedType(query)
            query._name = 'ns1:Tag'
            matrices = jaspar.searchByTag(query, format, database)
            for matrix in matrices:
                _print_matrix(matrix)

    except Exception, e:
            print 'Exception:', e            

    
def _print_matrix(matrix):
    
    print "Matrix id: " + matrix._attrs[(None, 'ID')]
    Aline = "A "
    for x in matrix.A.col:
        Aline += x + "\t" 
    print Aline
    Tline = "T "
    for x in matrix.T.col:
        Tline += x + "\t" 
    print Tline
    Cline = "C "
    for x in matrix.C.col:
        Cline += x + "\t" 
    print Cline
    Gline = "G "
    for x in matrix.G.col:
        Gline += x + "\t" 
    print Gline

    print "\nTags:"
    for x in matrix.Tag:
        if type(x.Value) == type([]):
            for v in x.Value:
                print x.Name + ": " + v
        else:
	    print x.Name + ": " + x.Value;
    print


if __name__ == '__main__':
 main ()



