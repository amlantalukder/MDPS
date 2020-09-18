import pdb, math, datetime
from utils import *

DIAGONAL = 0
VERTICAL = 1
HORIZONTAL = 2

mode = 1

INF = -float('inf')
pwm_matrix = []
trans_matrix = []
posTransMatrix = []

costDict = {}

output_data = ""

mrna_seq = ""
mirna_seq = ""
mrna_name = ""
mirna_name = ""
bond_str = ""
msg = ""

OUTPUT_FOLDER = ""

global_alignments = []

#-----------------------------------------------
def matchMaps(map1, map2):

    match_counter = 0

    total_match_pos = 0
    
    for item in map1:

        if item[1] == -1:
            continue

        found = False

        for item1 in map2:

            if item1[1] == -1:
                continue

            if item == item1:

                found = True
                break

        if found:
            match_counter += 1

        total_match_pos += 1

    return float(match_counter * 100/total_match_pos)

#-----------------------------------------------
def writeDataWithMarkerAsHTMLTable(data, filename):

    html = ""

    for items in data:
        tds = ""
        for item in items:

            imgName = ""

            if item[-2:] == "*\\":

                imgName = "red_angle_"
                item = item[:-2]

            elif item[-2:] == "*-":

                imgName = "red_left_"
                item = item[:-2]

            elif item[-2:] == "*|":

                imgName = "red_up_"
                item = item[:-2]


            if item[-1:] == "\\":

                imgName = "blue_angle_" + imgName
                item = item[:-1]

            elif item[-1:] == "-":

                imgName = "blue_left_" + imgName
                item = item[:-1]

            elif item[-1:] == "|":

                imgName = "blue_up_" + imgName
                item = item[:-1]


            if imgName != "":
                imgName = "<img src='example_files/" + imgName[:-1] + ".jpg" + "' height='20px'><br />"
            
            tds += "<td>" + imgName + str(item) + "</td>"
            
        html += "<tr>" + tds + "</tr>"

    html = "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 3.2//EN\"><html><head><script src=\"debug.js\" type=\"text/javascript\"></script></head>" + \
           "<body><div><img src=\"example_files/red_left.jpg\" height=\"20px\">&nbsp;CLASH alignment path<br><img src=\"example_files/blue_left.jpg\" height=\"20px\">&nbsp;DP alignment path</div><br>" + \
           "Please click on each cell to see the 3rd dimension of bulges<br /><br /><table style=\"width:100%\" border=\"1\">" + html + "</table>" + \
           "<div id=\"container\" style=\"width:500px; height:450px; padding:5px; z-index:3000; visibility:hidden; display:none; position:fixed; background-color:#ffffff; border:1px solid black;\">" + \
           "<input type=\"button\" style=\"float:right;\" value=\"close\" onclick=\"hideDiv(\'container\');\" /><div style=\"clear:both\"></div>" + \
           "<div id=\"containerInner\" style=\"width=100%; height:400px; overflow-y:scroll; border: 1px solid black; padding:10\"></div></div></body></html>"

    writeFile("../visualization/temp/{}.html".format(filename), html, "w")

#-----------------------------------------------
def writeDirectionArray(alignment, alignment_true=[], file_name='alignment'):

    global debug, seq1, seq2

    path = getAlignmentPath(getAlignmentDP())

    if len(alignment_true) > 0:
        path1 = getAlignmentPath(alignment_true)

    out = [[""] + [item for item in seq1]]
    out1 = [[""] + [item for item in seq1]]
    debug_str = ""

    for i in range(0, len(debug)):

        row = [seq2[i]]
        row1 = [seq2[i]]

        for j in range(0, len(debug[i])):
            
            pop_up_str = "<strong>Current State: s(" + str(i) + "," + str(j) + ")</strong><br /><br />Best alignment of:<br />" + seq1[j:0:-1] + "<br />" + seq2[i:0:-1] + \
                       "<br /><br />Costs:<br />Stack: " + str(stack[i][j][0]) + ", prev state: " + str(stack[i][j][1]) + "<br />Bulge (miRNA): " + str(bulge_mirna[i][j][0]) + \
                       ", prev state: " + str(bulge_mirna[i][j][1]) + "<br />Bulge (mRNA): " + str(bulge_mrna[i][j][0]) + ", prev state: " + str(bulge_mrna[i][j][1])

            dp_value = "<a onclick=\"showText(\'" + pop_up_str + "\')\">[" + str(i) + ", " + str(j) + "]</a>"
            
            if (i, j) in path:
                item = path[(i, j)]
                
                if len(alignment_true) > 0 and (i, j) in path1:
                    item1 = path1[(i, j)]
                    dp_value += ("<br />[" + item1[1] + "]<br />" + item[0] + "*" + item1[0])
                else:
                    dp_value += ("<br />[" + item[1] + "]<br />" + item[0])

            elif len(alignment_true) > 0 and (i, j) in path1:
                item1 = path1[(i, j)]
                dp_value += ("<br />[" + item1[1] + "]<br />*" + item1[0])

            row.append(dp_value)

            #dp_value = "S({}, {}, 0): ".format(i, j) + str(round(stack[i][j][0], 2)) + ", prev state: " + str(stack[i][j][1]) + \
             #          ", S({}, {}, 1): ".format(i, j) + str(round(bulge_mirna[i][j][0], 2)) + ", prev state: " + str(bulge_mirna[i][j][1]) + \
              #         ", S({}, {}, 2): ".format(i, j) + str(round(bulge_mrna[i][j][0], 2)) + ", prev state: " + str(bulge_mrna[i][j][1])
            dp_value = round(max([stack[i][j][0], bulge_mirna[i][j][0], bulge_mrna[i][j][0]]), 2)
            row1.append(dp_value)

        out.append(row)
        out1.append(row1)

    print(formatDataTable(out1))    

    out[0][0] = "&nbsp;"
    out[0][1] = "&nbsp;"
    out[1][0] = "&nbsp;"
        
    writeDataWithMarkerAsHTMLTable(out, file_name)

#-----------------------------------------------
def getStandardAlignmentCost(alignment):

    len_mirna = len([item for item in alignment if item[1] >= 0])

    cost = 0
    record = False

    for item in alignment[::-1]:

        if item[1] != -1:
            record = True

        if record:

            #pdb.set_trace()
            
            bond_type = item[2]

            if bond_type == BULGE_MRNA:
                cost += trans(prev_bond_type, bond_type)
                prev_bond_type = item[2]
            else:
                if item[1] == 0:
                    cost += math.log(pwm_matrix[bond_type][item[1]])
                else:
                    cost += trans(prev_bond_type, bond_type) + math.log(pwm_matrix[bond_type][item[1]])

                prev_bond_type = bond_type

        if item[1] == len_mirna-1:
            break

    return cost/len_mirna

#-----------------------------------------------
def getAlignments(arr):

    for item in arr:

        str1, str2, str3 = getAlignmentDP()

        print(str(item[1]) + "\n")

        print(str1 + "\n" + str3 + "\n" + str2)

#------------------------------------------------
def getAlignmentCostAndStartIndex():

    max_row = i = len(seq2)-1
    max_col = len(seq1)-1
    max_val = INF
    prev_direction = -1

    for j in range(len(seq1)-1, len(seq1)-2, -1):

        if stack[i][j][0] > max_val:

            max_val = stack[i][j][0]
            max_col = j
            prev_direction = 0

        if bulge_mirna[i][j][0] > max_val:

            max_val = bulge_mirna[i][j][0]
            max_col = j
            prev_direction = 1

        if bulge_mrna[i][j][0] > max_val:

            max_val = bulge_mrna[i][j][0]
            max_col = j
            prev_direction = 2
    
    return max_val/(len(seq2)-1), (max_row, max_col, prev_direction)

#-----------------------------------------------
def getAlignmentDP():
   
    max_row, max_col, prev_direction = getAlignmentStartIndex()

    alignment = []

    j = len(seq1)-1
    while j > max_col:
        alignment.append([len(mrna_seq)-j, -1, BULGE_MRNA])
        j -= 1

    i = max_row

    while(True):

        if prev_direction == -1:
            break
 
        if prev_direction == DIAGONAL:
    
            alignment.append([len(mrna_seq)-j, i-1, getBondType(seq1[j], seq2[i])])

            prev_direction = stack[i][j][1]

            i -= 1
            j -= 1

        elif prev_direction == VERTICAL:

            alignment.append([-1, i-1, BULGE_MIRNA])
            prev_direction = bulge_mirna[i][j][1]
            
            i -= 1
                        
        elif prev_direction == HORIZONTAL:

            alignment.append([len(mrna_seq)-j, -1, BULGE_MRNA])
            prev_direction = bulge_mrna[i][j][1]
            j -= 1
            

        else:
            pdb.set_trace()

    alignment += [[len(mrna_seq)-k, -1, BULGE_MRNA] for k in range(j, 0, -1)]

    return alignment

#-----------------------------------------------
def getAlignmentCost():

    cost, index = getAlignmentCostAndStartIndex()

    return cost

#-----------------------------------------------
def getAlignmentStartIndex():

    cost, index = getAlignmentCostAndStartIndex()

    return index

#-----------------------------------------------
def getCostMsg(pos_info, prev_state, total_cost):

    i, j, state = pos_info

    if state == BULGE_MRNA:

        if j > 0 and j < len(seq2) - 1:
            return str(trans(prev_state, BULGE_MRNA)) + " = " + str(total_cost)
        else:
            return str(total_cost)

    else:

        if j == 0:
            return str(pwm_matrix[state][j]) + " = " + str(total_cost)
        else:
            return str(trans(prev_state, state)) + " x " + str(pwm_matrix[state][j]) + " = " + str(total_cost)

    return ''
        
#-----------------------------------------------
def getAlignmentPath(alignment):

    path = {}

    alignment = alignment[::-1]

    record = False

    i, j = len(seq2)-1, len(seq1)-1
   
    for index in range(len(alignment)):

        item = alignment[index]

        if item[1] != -1:
            record = True

        if record:

            state = item[2]
            if index > 0:
                prev_state = alignment[index-1][2]
            else:
                prev_state = -1

            if state in [MATCH, WOBBLE, MISMATCH]:

                i, j = item[1]+1, len(seq1)-item[0]-1
                msg = getCostMsg(item, prev_state, stack[i][j][0])

                path[(i, j)] = ["\\", msg]

            elif state == BULGE_MIRNA:

                i = item[1]+1
                msg = getCostMsg(item, prev_state, bulge_mirna[i][j][0])

                path[(i, j)] = ["|", msg]
                
            elif state == BULGE_MRNA:

                j = len(seq1)-item[0]-1
                msg = getCostMsg(item, prev_state, bulge_mrna[i][j][0])

                path[(i, j)] = ["-", msg]

            else:
                pdb.set_trace()

        if item[1] == len(mirna_seq):
            break

    return path


#-----------------------------------------------
def loadDPArray():

    global stack, bulge_mirna, bulge_mrna, debug

    stack = [[[INF, -1] for j in range(0, len(seq1))] for i in range(0, len(seq2))]
    bulge_mirna = [[[INF, -1] for j in range(0, len(seq1))] for i in range(0, len(seq2))]
    bulge_mrna = [[[INF, -1] for j in range(0, len(seq1))] for i in range(0, len(seq2))]
    
    debug = [["" for j in range(0, len(seq1))] for i in range(0, len(seq2))]

    seq1_length = len(seq1)
    seq2_length = len(seq2)
    
    for i in range(0, seq2_length):
        for j in range(0, seq1_length):
            
            if i > 0 and j > 0: setCost(i, j, seq1[j], seq2[i])
            if i > 0: setCost(i, j, "-", seq2[i])
            if j > 1: setCost(i, j, seq1[j], "-")

#-----------------------------------------------
def setCost(i, j, x, y):

    if x == "-":

        state = BULGE_MIRNA

        #-----------------------------------------------
        # For first position of the miRNA, no transition
        # from previous state is needed
        #-----------------------------------------------
        if i == 1:
            prev_cost = 0
            prev_direction = -1
        #-----------------------------------------------
        # For all other positions, consider a transition
        # from previous bulge_mirna or stack state
        #-----------------------------------------------
        else:
            prev_stack_to_bulge_mirna = stack[i-1][j][0] + trans(getBondType(seq1[j], seq2[i-1]), state)
            prev_bulge_mirna_to_bulge_mirna = bulge_mirna[i-1][j][0] + trans(BULGE_MIRNA, state)

            if prev_stack_to_bulge_mirna > prev_bulge_mirna_to_bulge_mirna:
                prev_cost = prev_stack_to_bulge_mirna
                prev_direction = DIAGONAL
            else:
                prev_cost = prev_bulge_mirna_to_bulge_mirna
                prev_direction = VERTICAL
	#-----------------------------------------------
        # Add current state cost
        #-----------------------------------------------
        bulge_mirna[i][j][0] = prev_cost + math.log(pwm_matrix[state][i-1])
        bulge_mirna[i][j][1] = prev_direction

    elif y == "-":

        state = BULGE_MRNA

        #-----------------------------------------------
        # Consider a transition from previous
        # bulge_target or stack state
        #-----------------------------------------------
        prev_stack_to_bulge_mrna = stack[i][j-1][0] + trans(getBondType(seq1[j-1], seq2[i]), state)
        prev_bulge_mrna_to_bulge_mrna = bulge_mrna[i][j-1][0] + trans(BULGE_MRNA, state)

        if prev_stack_to_bulge_mrna > prev_bulge_mrna_to_bulge_mrna:
            prev_cost = prev_stack_to_bulge_mrna
            prev_direction = DIAGONAL
        else:
            prev_cost = prev_bulge_mrna_to_bulge_mrna
            prev_direction = HORIZONTAL

        bulge_mrna[i][j][0] = prev_cost
        bulge_mrna[i][j][1] = prev_direction

    else:

        state = getBondType(x, y)

        #-----------------------------------------------
        # For first position of the miRNA, no transition
        # from previous state is needed
        #-----------------------------------------------
        if i == 1:
            prev_cost = 0
            prev_direction = -1
        #-----------------------------------------------
        # For first position of the target, consider
        # transition from the previous bulge_mirna state
        #-----------------------------------------------
        elif j == 1:
            prev_cost = bulge_mirna[i-1][j-1][0] + trans(BULGE_MIRNA, state)
            prev_direction = VERTICAL
        #-----------------------------------------------
        # For all other positions, consider transition
        # from any of the previous states that costs least
        #-----------------------------------------------
        else:
            prev_stack_to_stack = stack[i-1][j-1][0] + trans(getBondType(seq1[j-1], seq2[i-1]), state)
            prev_bulge_mirna_to_stack = bulge_mirna[i-1][j-1][0] + trans(BULGE_MIRNA, state)
            prev_bulge_mrna_to_stack = bulge_mrna[i-1][j-1][0] + trans(BULGE_MRNA, state)

            if prev_stack_to_stack >= prev_bulge_mirna_to_stack:
                if prev_stack_to_stack >= prev_bulge_mrna_to_stack:
                    prev_cost = prev_stack_to_stack
                    prev_direction = DIAGONAL
                else:
                    prev_cost = prev_bulge_mrna_to_stack
                    prev_direction = HORIZONTAL
            else:
                if prev_bulge_mirna_to_stack >= prev_bulge_mrna_to_stack:
                    prev_cost = prev_bulge_mirna_to_stack
                    prev_direction = VERTICAL
                else:
                    prev_cost = prev_bulge_mrna_to_stack
                    prev_direction = HORIZONTAL

        #-----------------------------------------------
        # Add current state cost
        #-----------------------------------------------
        stack[i][j][0] = prev_cost + math.log(pwm_matrix[state][i-1])
        stack[i][j][1] = prev_direction

#-----------------------------------------------
def trans(prev_state, currState):

    transCost = math.log(trans_matrix[prev_state][currState])

    if math.isnan(transCost):

        return INF

    return transCost

#-----------------------------------------------
def matched(x, y):

    if (x, y) in [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("G", "T"), ("T", "G")]:
        return True

    return False

#-----------------------------------------------
def getBondType(x, y):

    if (x, y) in [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")]:
        return MATCH

    if (x, y) in [("G", "T"), ("T", "G")]:
        return WOBBLE
        #return MATCH

    return MISMATCH

#-----------------------------------------------
def sequenceAlignment():

    global seq1, seq2

    seq1 = "-" + mrna_seq[::-1]
    seq2 = "-" + mirna_seq

    #msg.append(["loadDPArray", datetime.datetime.now()])

    loadDPArray()

#-----------------------------------------------
def compareAlignments(mrna_seq, mirna_seq, alignment_true):

    global seq1, seq2

    seq1 = "-" + mrna_seq[::-1]
    seq2 = "-" + mirna_seq

    loadDPArray()

    alignment = getAlignmentDP()
    
    return matchMaps(alignment_true, alignment)

#-----------------------------------------------
def sequenceAlignmentCost():

    global seq1, seq2

    seq1 = "-" + mrna_seq[::-1]
    seq2 = "-" + mirna_seq[::-1]

    loadDPArray()

    return getAlignmentCost()
