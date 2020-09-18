import dynamic_programming_seq_alignment as dp
from utils import *

# -------------------------------------------------------------------------
def extendTargetPos(target_seq, start, end, size=40):

    if end-start+1 >= size:
        return start, end

    if start <= size/2:
        start = 1
    else:
        start -= size/2

    ext = size - (end-start+1)

    if end + ext > len(target_seq):
        end = len(target_seq)
    else:
        end += ext

    return start, end

# -------------------------------------------------------------------------
def getPredictionDict(prediction_data):

    prediction_dict = {}
    
    for item in prediction_data:

        mirna_id, target_id, start, end = item[0], item[1], int(item[2]), int(item[3])
        other_info = []
        if len(item) > 4: other_info = item[4:]
        start, end = extendTargetPos(target_db[target_id], start, end)

        if (target_id, mirna_id) not in prediction_dict:
            prediction_dict[(target_id, mirna_id)] = [[start, end] + other_info]
        else:
            prediction_dict[(target_id, mirna_id)].append([start, end] + other_info)

    return prediction_dict

# -------------------------------------------------------------------------
def checkMDPS(target_id, mirna_id, target_seq, mirna_seq):

    dp.mirna_name = mirna_id
    dp.mirna_seq = mirna_seq
    dp.mrna_name = target_id
    dp.mrna_seq = target_seq

    cost = dp.sequenceAlignmentCost()

    return cost >= (avg_cost - 2 * std_dev_cost)

# -------------------------------------------------------------------------
def predict(prediction_dict):

    filtered_predictions = []
    
    # ---------------------------------------------------------------------
    # Scan each predicted targets with MDPS
    # ---------------------------------------------------------------------
    print "Scanning the predicted positives..."
    perc = 10
    counter = 0

    count_original_predictions = 0
    
    for (target_id, mirna_id) in prediction_dict:

        counter += 1
        perc = showPercBar(counter, len(prediction_dict), perc)

        mirna_seq = mirna_db[mirna_id]

        if len(mirna_seq) > len(dp.pwm_matrix[0]):
            continue

        sites = prediction_dict[(target_id, mirna_id)]

        for info in sites:

            count_original_predictions += 1

            start, end = info[:2]

            target_seq = target_db[target_id][start-1:end]
                        
            if checkMDPS(target_id, mirna_id, target_seq, mirna_seq):
                filtered_predictions.append([mirna_id, target_id] + info)

    print '{}% of the original predictions was filtered by MDPS...'.format(100-len(filtered_predictions)*100./count_original_predictions)
    writeDataTableAsText(filtered_predictions, 'predictions_filtered.txt')

# -------------------------------------------------------------------------
def showArgError(parser):
    parser.parse_args(['-h'])
    exit()

# -------------------------------------------------------------------------
def showError(msg):
    print('----------- Error !!! ------------')
    print(msg)
    print('----------------------------------')
        
parser = argparse.ArgumentParser(description='Filter predicted miRNA targets', \
                                 epilog='Example: python mdps.py -p examples/test_predictions.txt -t examples/test_target_seq.fa -m examples/test_mirna_seq.fa -o examples/test_output.txt')
required_arg = parser.add_argument_group('required arguments')
required_arg.add_argument('-p', help="Path for miRNA prediction file in a tsv format", required=True, metavar="PREDICTIONS")
required_arg.add_argument('-t', help='Path for target sequences in fasta file format', required=True, metavar="TARGET_SEQUENCES")
required_arg.add_argument('-m', help='Path for miRNA sequences in fasta file format', required=True, metavar="MIRNA_SEQUENCES")
parser.add_argument('-o', help='Path for MDPS outputs', metavar="OUTPUT")

if '-h' in sys.argv[1:]:
    showArgError(parser)

opts,args=getopt.getopt(sys.argv[1:],'p:t:m:o:')

prediction_file_path = ''
target_seq_path = ''
mirna_seq_path = ''
output_path = 'predictions_filtered.txt'

for i in opts:
    if i[0] == '-p':
        prediction_file_path = i[1]
    elif i[0] == '-t':
        target_seq_path = i[1]
    elif i[0] == '-m':
        mirna_seq_path = i[1]
    elif i[0] == '-o':
        output_path = i[1]
    else:
        showArgError(parser)

if prediction_file_path == '' or target_seq_path == '' or mirna_seq_path == '' or output_path == '':
    showArgError(parser)
    
# Validate gene path
if not os.path.isfile(prediction_file_path):
    showError('The following file path for miRNA target predictions does not exist' + '\n' + prediction_file_path)
    exit()

# Validate chromosome sizes file path
if not os.path.isfile(target_seq_path):
    showError('The following file path for target sequences does not exist' + '\n' + target_seq_path)
    exit()

# Validate chromosome sizes file path
if not os.path.isfile(mirna_seq_path):
    showError('The following file path for miRNA sequences does not exist' + '\n' + mirna_seq_path)
    exit()


print 'Getting MDPS model...'

position_wise_knowledge_dict = ast.literal_eval(readFile('parameters/mirna_position_wise_knowledge.txt')[0])
dp.pwm_matrix = [[float(item1) if item1 != 'None' else None for item1 in item] for item in position_wise_knowledge_dict['common'][0]]
dp.trans_matrix = [[float(item1) if item1 != 'None' else None for item1 in item] for item in position_wise_knowledge_dict['common'][1]]

average_cost_dict = ast.literal_eval(readFile('parameters/cost_threshold.txt')[0])
cost_info = average_cost_dict['common']
avg_cost, max_cost, min_cost, std_dev_cost = float(cost_info[0])/float(cost_info[1]), \
                                             float(cost_info[2]), float(cost_info[3]), \
                                             float(cost_info[4])

print 'Reading target sequences...'
target_db = fastaToDict(target_seq_path)

print 'Reading miRNA sequences...'
mirna_db = fastaToDict(mirna_seq_path)

print 'Reading predicted data...'
ext_predictions = readFileInTable(prediction_file_path)[:100]
ext_prediction_dict = getPredictionDict(ext_predictions)

print 'Scoring predictions...'
predict(ext_prediction_dict)

print 'Done...'
