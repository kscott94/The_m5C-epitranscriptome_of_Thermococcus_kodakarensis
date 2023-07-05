def load_genome(path):
    with open(path, 'r') as genome:
        """remove fasta header"""
        genome_sequence_initial = genome.readlines()[1:]
        """genome sequence as string. Add an "X" to the front to correct for zero=based indexing"""
        genome_sequence = "X" + ''.join(genome_sequence_initial).replace('\n', '')
        return(genome_sequence)


def rev_comp(seq, type, convert = True):
    #seq; nucleotide sequence
    #type; DNA or RNA to convert to

    seq = seq.upper()
    if convert == True and type == 'RNA':
        seq = seq.replace("T", "U")
    if convert == True and type == 'DNA':
        seq = seq.replace("U", "T")

    dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    if type=='RNA':
        return ''.join([rna_complement[base] for base in seq[::-1]])
    elif type == 'DNA':
        return ''.join([dna_complement[base] for base in seq[::-1]])


def convert(seq, type):
    #convert DNA to RNA or vis versa
    seq = seq.upper()
    if type == 'RNA':
        seq = seq.replace("T", "U")
    if type == 'DNA':
        seq = seq.replace("U", "T")
    return(seq)


def ref_strand(nuc):
    """
    determine if the strand is + ir -
    :param nuc: C or G
    :return: + or -
    """
    #nuc: C or G
    nuc=nuc.upper()
    if nuc=='G':
        ref = '-'
    elif nuc=='C':
        ref = '+'
    else:
        ref = 'NA'
    return(ref)


def get_descript(description):
    """
    make a dictionary of tags in annotation file where the attribute is the key and the value is the... value
    :param description: description line in annotation file.
    :return: dictionary of key/value pairs of tags in annotation description
    """
    desc_dict = {}
    ls = description.split(";")
    for i in ls:
        name, value = i.split("=")
        desc_dict[name]=value
    return(desc_dict)


def codon_pos(pos):
    """
    get the codon position (1,2, or 3) from the position in transcript
    :param pos: position from start of transcript
    :return: codon position
    """
    if pos % 3 == 1:
        return 1
    if pos % 3 == 2:
        return 2
    if pos % 3 == 0:
        return 3
    else:
        return 'NaN'


def codon_seq(codon_pos, surroudning_seq):
    #codon_pos is 1,2, or 3
    #surrouding_seq is the sequence 2nt up and down stream

    aa_seq = str(surroudning_seq).replace("T", "U")
    if codon_pos == 1:
        return aa_seq[2:]
    elif codon_pos == 2:
        return aa_seq[1:4]
    elif codon_pos == 3:
        return aa_seq[:3]


def aa_id(aa_seq):
    """
    get the amino acid one letter symbol (ID) from the 3-letter codon sequence
    :param aa_seq: 3 letter codon sequence
    :return: amino acid ID
    """
    aa_seq = aa_seq.upper().replace("T","U")
    genecode = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W'}
    return(genecode[aa_seq])


def get_aaseq(seq):
    """
    convert a gene sequence into an amino acid sequence
    :param seq: gene sequence (CDS)
    :return: list of amino acids that make up gene
    """
    start = 0
    stop = 3
    iterations = int(len(seq)/3) - 1
    aa_sequence = ''
    for i in range(iterations):
        aaseq = seq[start:stop]
        aaid = aa_id(aaseq)
        aa_sequence += aaid
        start +=3
        stop += 3
    aa_sequence = aa_sequence.replace("T", "U")
    return(aa_sequence)

def get_aacodon(seq):
    """
    convert a gene sequence into a codon sequence
    :param seq: gene sequence (CDS)
    :return: list of codons that make up gene
    """
    start = 0
    stop = 3
    iterations = int(len(seq)/3) - 1
    aa_sequence = []
    for i in range(iterations):
        aaseq = seq[start:stop]
        aa_sequence.append(aaseq)
        start +=3
        stop += 3
    return(aa_sequence)


def ls_to_csv(ls):
    if len(ls) == 0:
        return("NA")
    string = None
    for i in ls:
        if string is None:
            string = str(i) + ","
        else:
            string = string + str(i) + ","

    string = string.rstrip(",")
    return(string)


def get_KOD1_position(query, subject):
    search_length = len(query)
    start_search = 1
    stop_search = start_search + search_length
    all_positions = []

    while stop_search <= len(subject):
        if query == subject[start_search:stop_search]:
            all_positions.append(start_search)

        start_search += 1
        stop_search = start_search + search_length

    if len(all_positions) == 0:
        return('NA')
    else:
        pos = ''.join(f'{i};'for i in all_positions).rstrip(';')
        return(pos)


def convert_tag(tag, output='old_locus_tag'):
    dict = {
        "TK_RS00370": "TK_RS00370", "TK_RS02665": "TK_RS02665", "TK_RS11710": "TK_RS11710", "TK_RS11730": "TK_RS11730",
        "TK_RS00080": "TK0015", "TK_RS00115": "TK0022", "TK_RS00145": "TK0028", "TK_RS00155": "TK0030",
        "TK_RS00180": "TK0034", "TK_RS00185": "TK0035", "TK_RS00190": "TK0036", "TK_RS00195": "TK0037",
        "TK_RS00200": "TK0038", "TK_RS00210": "TK0040", "TK_RS00220": "TK0042", "TK_RS00230": "TK0044",
        "TK_RS00240": "TK0046", "TK_RS00290": "TK0057", "TK_RS00300": "TK0059", "TK_RS00355": "TK0070",
        "TK_RS00530": "TK0107", "TK_RS00555": "TK0112", "TK_RS00580": "TK0117", "TK_RS00615": "TK0124",
        "TK_RS00625": "TK0126", "TK_RS00635": "TK0129", "TK_RS00640": "TK0130", "TK_RS00645": "TK0131",
        "TK_RS00650": "TK0132", "TK_RS00665": "TK0135", "TK_RS00695": "TK0142", "TK_RS00700": "TK0143",
        "TK_RS00735": "TK0150", "TK_RS00760": "TK0155", "TK_RS00765": "TK0156", "TK_RS00790": "TK0161",
        "TK_RS00800": "TK0163", "TK_RS00805": "TK0164", "TK_RS00815": "TK0166", "TK_RS00895": "TK0182",
        "TK_RS00900": "TK0183", "TK_RS00915": "TK0186", "TK_RS00930": "TK0189", "TK_RS00955": "TK0194",
        "TK_RS00960": "TK0195", "TK_RS01110": "TK0226", "TK_RS01170": "TK0238", "TK_RS01390": "TK0284",
        "TK_RS01450": "TK0298", "TK_RS01460": "TK0299", "TK_RS01485": "TK0304", "TK_RS01505": "TK0306",
        "TK_RS01515": "TK0307", "TK_RS01520": "TK0308", "TK_RS01525": "TK0309", "TK_RS01545": "TK0313",
        "TK_RS01615": "TK0327", "TK_RS01630": "TK0330", "TK_RS01715": "TK0348", "TK_RS01725": "TK0350",
        "TK_RS01735": "TK0352", "TK_RS01865": "TK0378", "TK_RS01975": "TK0400", "TK_RS02085": "TK0423",
        "TK_RS02100": "TK0426", "TK_RS02105": "TK0427", "TK_RS02175": "TK0441", "TK_RS02180": "TK0442",
        "TK_RS02190": "TK0444", "TK_RS02200": "TK0446", "TK_RS02235": "TK0453", "TK_RS02275": "TK0461",
        "TK_RS02295": "TK0465", "TK_RS02305": "TK0467", "TK_RS02325": "TK0471", "TK_RS02345": "TK0475",
        "TK_RS02370": "TK0480", "TK_RS02375": "TK0481", "TK_RS02435": "TK0493", "TK_RS02440": "TK0494",
        "TK_RS02455": "TK0495", "TK_RS02470": "TK0499", "TK_RS02515": "TK0508", "TK_RS02525": "TK0510",
        "TK_RS02595": "TK0524", "TK_RS02600": "TK0525", "TK_RS02615": "TK0528", "TK_RS02635": "TK0535",
        "TK_RS02640": "TK0536", "TK_RS02645": "TK0537", "TK_RS02670": "TK0541", "TK_RS02700": "TK0548",
        "TK_RS02710": "TK0550", "TK_RS02740": "TK0556", "TK_RS02750": "TK0558", "TK_RS02760": "TK0560",
        "TK_RS02805": "TK0569", "TK_RS02810": "TK0570", "TK_RS02820": "TK0572", "TK_RS02830": "TK0574",
        "TK_RS02915": "TK0589", "TK_RS03040": "TK0615", "TK_RS03045": "TK0616", "TK_RS03110": "TK0629",
        "TK_RS03175": "TK0643", "TK_RS03195": "TK0647", "TK_RS03205": "TK0649", "TK_RS03210": "TK0650",
        "TK_RS03230": "TK0654", "TK_RS03240": "TK0656", "TK_RS03245": "TK0657", "TK_RS03270": "TK0662",
        "TK_RS03275": "TK0663", "TK_RS03285": "TK0665", "TK_RS03300": "TK0668", "TK_RS03305": "TK0669",
        "TK_RS11630": "TK0675", "TK_RS11635": "TK0676", "TK_RS03350": "TK0678", "TK_RS03370": "TK0682",
        "TK_RS03375": "TK0683", "TK_RS03380": "TK0684", "TK_RS03485": "TK0705", "TK_RS03490": "TK0706",
        "TK_RS03510": "TK0710", "TK_RS03640": "TK0736", "TK_RS03710": "TK0751", "TK_RS03750": "TK0759",
        "TK_RS03785": "TK0765", "TK_RS03820": "TK0773", "TK_RS03860": "TK0781", "TK_RS03925": "TK0793",
        "TK_RS04010": "TK0810", "TK_RS04015": "TK0811", "TK_RS04030": "TK0814", "TK_RS04040": "TK0816",
        "TK_RS04085": "TK0826", "TK_RS04105": "TK0830", "TK_RS04170": "TK0844", "TK_RS04195": "TK0850",
        "TK_RS04200": "TK0853", "TK_RS04265": "TK0865", "TK_RS04270": "TK0866", "TK_RS04295": "TK0871",
        "TK_RS04330": "TK0878", "TK_RS04365": "TK0885", "TK_RS04410": "TK0894", "TK_RS04415": "TK0895",
        "TK_RS04420": "TK0896", "TK_RS04440": "TK0900", "TK_RS04445": "TK0901", "TK_RS04450": "TK0902",
        "TK_RS04455": "TK0903", "TK_RS04495": "TK0911", "TK_RS04510": "TK0914", "TK_RS04520": "TK0916",
        "TK_RS04530": "TK0918", "TK_RS04560": "TK0925", "TK_RS04605": "TK0935", "TK_RS04635": "TK0942",
        "TK_RS04640": "TK0943", "TK_RS04645": "TK0944", "TK_RS04655": "TK0946", "TK_RS04700": "TK0955",
        "TK_RS04720": "TK0959", "TK_RS04755": "TK0966", "TK_RS04760": "TK0967", "TK_RS04765": "TK0968",
        "TK_RS04800": "TK0975", "TK_RS04835": "TK0982", "TK_RS04850": "TK0985", "TK_RS04855": "TK0986",
        "TK_RS04880": "TK0991", "TK_RS04935": "TK1002", "TK_RS04940": "TK1003", "TK_RS04945": "TK1004",
        "TK_RS04950": "TK1005", "TK_RS04955": "TK1006", "TK_RS05000": "TK1016", "TK_RS05010": "TK1018",
        "TK_RS05025": "TK1021", "TK_RS05040": "TK1024", "TK_RS05045": "TK1025", "TK_RS05115": "TK1039",
        "TK_RS05120": "TK1040", "TK_RS11655": "TK1046", "TK_RS05160": "TK1048", "TK_RS05170": "TK1050",
        "TK_RS05175": "TK1051", "TK_RS05200": "TK1056", "TK_RS05215": "TK1060", "TK_RS05230": "TK1063",
        "TK_RS05290": "TK1077", "TK_RS05330": "TK1085", "TK_RS05345": "TK1088", "TK_RS05375": "TK1094",
        "TK_RS05400": "TK1099", "TK_RS05415": "TK1102", "TK_RS05440": "TK1107", "TK_RS05445": "TK1108",
        "TK_RS05450": "TK1109", "TK_RS05460": "TK1111", "TK_RS05465": "TK1112", "TK_RS05480": "TK1115",
        "TK_RS05510": "TK1120", "TK_RS05525": "TK1123", "TK_RS05575": "TK1133", "TK_RS05580": "TK1134",
        "TK_RS05620": "TK1140", "TK_RS05655": "TK1147", "TK_RS05710": "TK1158", "TK_RS05730": "TK1164",
        "TK_RS05745": "TK1167", "TK_RS05760": "TK1170", "TK_RS05770": "TK1172", "TK_RS05775": "TK1173",
        "TK_RS05795": "TK1177", "TK_RS05845": "TK1186", "TK_RS05870": "TK1191", "TK_RS05880": "TK1193",
        "TK_RS05965": "TK1209", "TK_RS06055": "TK1227", "TK_RS06065": "TK1229", "TK_RS06080": "TK1232",
        "TK_RS06085": "TK1233", "TK_RS06130": "TK1242", "TK_RS06145": "TK1245", "TK_RS06180": "TK1251",
        "TK_RS06195": "TK1254", "TK_RS06210": "TK1257", "TK_RS06220": "TK1259", "TK_RS06240": "TK1263",
        "TK_RS06300": "TK1274", "TK_RS06310": "TK1276", "TK_RS06340": "TK1282", "TK_RS06350": "TK1284",
        "TK_RS06355": "TK1285", "TK_RS06410": "TK1292", "TK_RS06425": "TK1295", "TK_RS06445": "TK1299",
        "TK_RS06470": "TK1304", "TK_RS06485": "TK1307", "TK_RS06490": "TK1308", "TK_RS06495": "TK1309",
        "TK_RS06505": "TK1311", "TK_RS06515": "TK1313", "TK_RS06530": "TK1316", "TK_RS06535": "TK1317",
        "TK_RS06550": "TK1320", "TK_RS06560": "TK1322", "TK_RS06585": "TK1327", "TK_RS06605": "TK1331",
        "TK_RS06680": "TK1344", "TK_RS06845": "TK1379", "TK_RS06900": "TK1390", "TK_RS06970": "TK1404",
        "TK_RS06975": "TK1405", "TK_RS06990": "TK1408", "TK_RS07000": "TK1410", "TK_RS07015": "TK1413",
        "TK_RS07030": "TK1415", "TK_RS07035": "TK1416", "TK_RS07045": "TK1418", "TK_RS07070": "TK1423",
        "TK_RS07110": "TK1428", "TK_RS07125": "TK1431", "TK_RS07235": "TK1454", "TK_RS07245": "TK1456",
        "TK_RS07265": "TK1460", "TK_RS07275": "TK1462", "TK_RS07280": "TK1463", "TK_RS07305": "TK1468",
        "TK_RS07315": "TK1470", "TK_RS07365": "TK1480", "TK_RS07380": "TK1481", "TK_RS07385": "TK1482",
        "TK_RS07395": "TK1484", "TK_RS07405": "TK1487", "TK_RS07415": "TK1489", "TK_RS07435": "TK1492",
        "TK_RS07445": "TK1494", "TK_RS07450": "TK1495", "TK_RS07485": "TK1500", "TK_RS07510": "TK1504",
        "TK_RS07555": "TK1513", "TK_RS07580": "TK1518", "TK_RS07585": "TK1519", "TK_RS07670": "TK1536",
        "TK_RS07720": "TK1546", "TK_RS07730": "TK1548", "TK_RS07770": "TK1556", "TK_RS07775": "TK1557",
        "TK_RS07800": "TK1562", "TK_RS07840": "TK1571", "TK_RS07870": "TK1577", "TK_RS07875": "TK1578",
        "TK_RS07895": "TK1582", "TK_RS07975": "TK1598", "TK_RS07980": "TK1599", "TK_RS08000": "TK1603",
        "TK_RS08005": "TK1604", "TK_RS08090": "TK1621", "TK_RS08105": "TK1624", "TK_RS08115": "TK1626",
        "TK_RS08150": "TK1633", "TK_RS08170": "TK1637", "TK_RS08185": "TK1640", "TK_RS08200": "TK1643",
        "TK_RS08225": "TK1648", "TK_RS08250": "TK1653", "TK_RS08275": "TK1658", "TK_RS08315": "TK1666",
        "TK_RS08375": "TK1679", "TK_RS08380": "TK1680", "TK_RS08400": "TK1684", "TK_RS08420": "TK1688",
        "TK_RS08445": "TK1691", "TK_RS08455": "TK1693", "TK_RS08460": "TK1694", "TK_RS08465": "TK1695",
        "TK_RS08470": "TK1696", "TK_RS08490": "TK1700", "TK_RS08565": "TK1713", "TK_RS08630": "TK1726",
        "TK_RS08680": "TK1736", "TK_RS08685": "TK1737", "TK_RS08690": "TK1738", "TK_RS08695": "TK1739",
        "TK_RS08705": "TK1741", "TK_RS08710": "TK1742", "TK_RS08720": "TK1744", "TK_RS08730": "TK1746",
        "TK_RS08740": "TK1748", "TK_RS08750": "TK1750", "TK_RS08765": "TK1753", "TK_RS08850": "TK1768",
        "TK_RS08855": "TK1769", "TK_RS08865": "TK1771", "TK_RS08880": "TK1774", "TK_RS08885": "TK1775",
        "TK_RS08895": "TK1777", "TK_RS08905": "TK1779", "TK_RS08915": "TK1781", "TK_RS08935": "TK1785",
        "TK_RS08940": "TK1786", "TK_RS08955": "TK1789", "TK_RS08975": "TK1793", "TK_RS08990": "TK1796",
        "TK_RS09010": "TK1800", "TK_RS09030": "TK1804", "TK_RS09100": "TK1818", "TK_RS09210": "TK1840",
        "TK_RS09220": "TK1842", "TK_RS09260": "TK1850", "TK_RS09285": "TK1855", "TK_RS09365": "TK1871",
        "TK_RS09405": "TK1879", "TK_RS09425": "TK1883", "TK_RS09445": "TK1887", "TK_RS09485": "TK1895",
        "TK_RS09505": "TK1899", "TK_RS09575": "TK1911", "TK_RS09580": "TK1912", "TK_RS09595": "TK1915",
        "TK_RS09660": "TK1929", "TK_RS09715": "TK1940", "TK_RS09745": "TK1946", "TK_RS09755": "TK1948",
        "TK_RS09760": "TK1949", "TK_RS09765": "TK1950", "TK_RS09770": "TK1951", "TK_RS09775": "TK1952",
        "TK_RS09780": "TK1953", "TK_RS09790": "TK1955", "TK_RS09805": "TK1958", "TK_RS09810": "TK1959",
        "TK_RS09860": "TK1966", "TK_RS09900": "TK1974", "TK_RS09920": "TK1978", "TK_RS09930": "TK1980",
        "TK_RS09935": "TK1981", "TK_RS09950": "TK1984", "TK_RS10025": "TK1999", "TK_RS10030": "TK2000",
        "TK_RS10035": "TK2001", "TK_RS10060": "TK2006", "TK_RS10130": "TK2021", "TK_RS10140": "TK2023",
        "TK_RS10145": "TK2024", "TK_RS10200": "TK2035", "TK_RS10210": "TK2037", "TK_RS10270": "TK2049",
        "TK_RS10285": "TK2052", "TK_RS10335": "TK2062", "TK_RS10340": "TK2063", "TK_RS10350": "TK2065",
        "TK_RS10370": "TK2069", "TK_RS10400": "TK2075", "TK_RS10410": "TK2077", "TK_RS10415": "TK2078",
        "TK_RS10420": "TK2079", "TK_RS10490": "TK2093", "TK_RS10525": "TK2100", "TK_RS10540": "TK2103",
        "TK_RS10545": "TK2104", "TK_RS10555": "TK2106", "TK_RS10580": "TK2110", "TK_RS10590": "TK2112",
        "TK_RS10600": "TK2114", "TK_RS10605": "TK2115", "TK_RS10615": "TK2117", "TK_RS10675": "TK2129",
        "TK_RS10695": "TK2131", "TK_RS10715": "TK2134", "TK_RS10725": "TK2136", "TK_RS10735": "TK2138",
        "TK_RS10745": "TK2140", "TK_RS10780": "TK2146", "TK_RS10855": "TK2160", "TK_RS10865": "TK2162",
        "TK_RS10875": "TK2164", "TK_RS10880": "TK2165", "TK_RS10935": "TK2176", "TK_RS10960": "TK2181",
        "TK_RS11030": "TK2195", "TK_RS11055": "TK2200", "TK_RS11075": "TK2204", "TK_RS11090": "TK2207",
        "TK_RS11140": "TK2217", "TK_RS11150": "TK2219", "TK_RS11175": "TK2224", "TK_RS11190": "TK2227",
        "TK_RS11215": "TK2232", "TK_RS11230": "TK2235", "TK_RS11255": "TK2237", "TK_RS11265": "TK2239",
        "TK_RS11290": "TK2244", "TK_RS11330": "TK2252", "TK_RS11350": "TK2256", "TK_RS11395": "TK2264",
        "TK_RS11425": "TK2268", "TK_RS11435": "TK2270", "TK_RS11440": "TK2271", "TK_RS11450": "TK2273",
        "TK_RS11460": "TK2275", "TK_RS11475": "TK2278", "TK_RS11515": "TK2286", "TK_RS11520": "TK2287",
        "TK_RS11530": "TK2289", "TK_RS11535": "TK2290", "TK_RS11545": "TK2292", "TK_RS11580": "TK2299",
        "TK_RS11600": "TK2303", "TK_RS00170": "TKr01", "TK_RS09830": "TKr04", "TK_RS11235": "TKr05",
        "TK_RS11245": "TKr06", "TK_RS00520": "TKt01", "TK_RS01495": "TKt02", "TK_RS01500": "TKt03",
        "TK_RS01510": "TKt04", "TK_RS01795": "TKt05", "TK_RS02445": "TKt07", "TK_RS02450": "TKt08",
        "TK_RS02835": "TKt09", "TK_RS03885": "TKt11", "TK_RS05600": "TKt12", "TK_RS05605": "TKt13",
        "TK_RS05830": "TKt14", "TK_RS06170": "TKt15", "TK_RS06370": "TKt16", "TK_RS06375": "TKt17",
        "TK_RS06380": "TKt18", "TK_RS07025": "TKt19", "TK_RS07085": "TKt20", "TK_RS07090": "TKt21",
        "TK_RS07095": "TKt22", "TK_RS07370": "TKt23", "TK_RS07430": "TKt25", "TK_RS07455": "TKt26",
        "TK_RS07475": "TKt27", "TK_RS07500": "TKt28", "TK_RS07525": "TKt29", "TK_RS08425": "TKt30",
        "TK_RS08430": "TKt31", "TK_RS08530": "TKt32", "TK_RS09545": "TKt34", "TK_RS09560": "TKt35",
        "TK_RS09835": "TKt36", "TK_RS09840": "TKt37", "TK_RS10710": "TKt39", "TK_RS10810": "TKt40",
        "TK_RS11240": "TKt42", "TK_RS11400": "TKt44", "TK_RS11405": "TKt45", "TK_RS06840": "TKt46"
    }

    if output == 'old_locus_tag' or 'TKID' in output.upper():
        out_tag = dict[tag]
        return(out_tag)

    elif output == 'locus_tag':
        for locus_tag, TKid in dict.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
            if tag == TKid:
                out_tag = locus_tag
                return (out_tag)
