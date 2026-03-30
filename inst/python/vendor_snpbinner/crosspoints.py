'''Uses genotyped SNP data to identify likely crossover points.
(See README for more information)'''
from math import log
import warnings
import os

def _crosspoints_batcher(input_path,output_path,predicted_homogeneity,predicted_cross_count,chrom_len,min_state_length=None,min_state_ratio=None):
    input_list = input_path
    if len(input_list) == 1:
        crosspoints(input_path[0],output_path,predicted_homogeneity,predicted_cross_count,chrom_len,min_state_length,min_state_ratio)
    else:
        output_folder = output_path
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for file in input_list:
            print(("Determining crosspoints for: "+file))
            out_file = os.path.join(output_folder,os.path.basename(file)+".crosp.csv")
            try:
                crosspoints(file,out_file,predicted_homogeneity,predicted_cross_count,0,min_state_length,min_state_ratio)
            except Exception as detail:
                print("FAILED! Details:")
                print(detail)

def crosspoints(input_path,output_path,predicted_homogeneity,predicted_cross_count,chrom_len,min_state_length=None,min_state_ratio=None):
    '''Runs the module'''

    #get file statistics
    individual_count,snp_count,auto_chrom_len = _get_file_stats(input_path)
    if chrom_len==0: 
        chrom_len=auto_chrom_len

    #determine the min distance between two crossover points
    if min_state_ratio!=None:
        min_state_length = chrom_len*float(min_state_ratio)

    # clear output file
    with open(output_path,"w") as outfile:
        pass

    # contruct the Hidden Markov Models (_HMM) used to predict states
    crosspoint_probability = predicted_cross_count/float(chrom_len)
    intra_p = predicted_homogeneity # emmision probability of the genotype corresponding to the state
    inter_p = 1-intra_p             # emmision probability of the opposite genotype

    # create a _HMM that registers heterogeneous regions
    error = crosspoint_probability/2.0  # transition probability for other 2 states
    crrct = 1-(error*2)                 # transition probability for same state
    
    # === FUZZY-HMM OVERRIDES ===
    try:
        _intra_override = os.environ.get('FSNP_INTRA_P', None)
        if _intra_override is not None:
            # 日志显示这里已经生效了！
            print(f"   [Core] Overriding intra_p: {intra_p:.4f} -> {float(_intra_override):.4f}") 
            intra_p = float(_intra_override)
            if intra_p < 1e-6: intra_p = 1e-6
            if intra_p > 1.0 - 1e-6: intra_p = 1.0 - 1e-6
            inter_p = 1.0 - intra_p

        _trans_mult = os.environ.get('FSNP_TRANS_MULT', None)
        if _trans_mult is not None:
            tm = float(_trans_mult)
            if tm < 0.05: tm = 0.05
            if tm > 20.0: tm = 20.0
            error = (crosspoint_probability/2.0) * tm
            if error > 0.49: error = 0.49
            crrct = 1.0 - (error * 2.0)
    except Exception as e:
        print(f"   [Core] Override warning: {e}")
    # === END OVERRIDES ===

    hmm_all = _HMM(
        states      = ["a", "b", "h"],
        priors      = [0.3, 0.3, 0.3],
        transition = [[crrct, error, error],
                      [error, crrct, error],
                      [error, error, crrct]],
        observable  = ["a", "b", "h"],
        emission   = [[intra_p, inter_p, inter_p],   
                      [inter_p, intra_p, inter_p],
                      [0.5,     0.5,     inter_p*2]] 
        )

    # create a _HMM that does not register heterogeneous regions
    error = crosspoint_probability
    crrct = 1-error
    hmm_nohet = _HMM(
        states      = ["a", "b"],
        priors      = [0.5, 0.5],
        transition = [[crrct, error],
                      [error, crrct]],
        observable  = ["a", "b"],
        emission   = [[intra_p, inter_p],
                      [inter_p, intra_p]] 
        )

    #Runs the crosspoint identification
    for i in range(0,individual_count):
        snplist,name = _read_column(input_path,i)
        if len(snplist) < 1:
            warnings.warn("HMM cannot be used on an empty data set, skippinng this line.", UserWarning, stacklevel=2)
            continue
        cp = _find_crosspoints(
            snplist          = snplist,
            min_state_length = min_state_length,
            chrom_length     = chrom_len,
            hmm_nohet         = hmm_nohet,
            hmm_all          = hmm_all)
        with open(output_path,"a") as outfile:
            outfile.write(",".join([name]+[str(n) for n in cp])+",\n")

def _find_crosspoints(snplist, min_state_length, chrom_length, hmm_nohet, hmm_all):
    '''Identifies crosspoints using two _HMMs.'''

    cross_points  = hmm_all.gapped_viterbi(snplist)
    no_hetero_cross_points = hmm_nohet.gapped_viterbi(snplist)

    #Replace any heterogenous regions which are below the minimum state length
    for i in range(2,len(cross_points)-1,2):
        if cross_points[i]-cross_points[i-2]<min_state_length:
            if cross_points[i-1] == hmm_all.states[-1] and (hmm_all.states[-1] not in (cross_points[i-3],cross_points[i+1])):
                start,stop = None,None
                for index,val in list(enumerate(no_hetero_cross_points))[0::2]:
                    if cross_points[i]<=val and stop==None:
                        stop = index+1
                        break
                    if cross_points[i-2]>=val:
                        start = index
                if not stop:
                    stop = no_hetero_cross_points[-1]
                cross_points[i-1:i] = no_hetero_cross_points[start+1:stop-1]

    #Identify all regions which are below the min state length
    too_short = []
    for i in range(2,len(cross_points),2):
        if cross_points[i]-cross_points[i-2]<min_state_length:
            if len(too_short)>0 and too_short[-1][-1][1] == i-2:
                too_short[-1].append((i-2,i))
            else:
                too_short.append([(i-2,i)])

    front_skipped,back_skipped = False,False
    if too_short and too_short[0][0][0]==0:
        front_skipped = True
        del too_short[0][0]
        if not too_short[0]:
            del too_short[0]
    if too_short and too_short[-1][-1][1]==len(cross_points)-1:
        back_skipped = True
        del too_short[-1][-1]
        if not too_short[-1]:
            del too_short[-1]

    for group in too_short:
        while len(group)>0:
            left_match = cross_points[group[0][0]-1] == cross_points[group[0][1]-1]
            right_match = cross_points[group[-1][1]+1] == cross_points[group[-1][1]-1] 

            if left_match or right_match: 
                if len(group)==1: 
                    break
                if left_match:  
                    group[:] = group[1:]
                if right_match: 
                    group[:] = group[:-1]
            else: 
                if cross_points[group[-1][1]]-cross_points[group[0][0]] >= min_state_length:
                    first_snp = None
                    last_snp = None
                    i = 0
                    while (first_snp == None or last_snp == None) and i<len(snplist):
                        if snplist[i][0] == cross_points[group[0][0]]: first_snp = i
                        if snplist[i][0] == cross_points[group[-1][1]]: last_snp = i
                        i+=1
                    group_snps = snplist[first_snp:last_snp+1]
                    most_likely = hmm_all.get_most_likely_state(group_snps)
                    for reg in group:
                        cross_points[reg[1]-1] = most_likely
                    group = []
                    continue
                
                if hmm_all.states[-1] in (cross_points[group[0][0]-1],cross_points[group[-1][1]+1]):
                    if cross_points[group[0][0]-1] == hmm_all.states[-1]: 
                        cross_points[group[0][1]-1] = hmm_all.states[-1]
                        group[:] = group[1:]
                        if len(group)==0: continue
                    if cross_points[group[-1][1]+1] == hmm_all.states[-1]: 
                        cross_points[group[-1][1]-1] = hmm_all.states[-1]
                        group[:] = group[:-1]
                else: 
                    left_size = cross_points[group[0][1]]-cross_points[group[0][0]]
                    right_size = cross_points[group[-1][1]]-cross_points[group[-1][0]]
                    if left_size<right_size:
                        cross_points[group[0][1]-1] = cross_points[group[0][0]-1]
                        group[:] = group[1:]
                    else:
                        cross_points[group[-1][1]-1] = cross_points[group[-1][1]+1]
                        group[:] = group[:-1]

    if len(cross_points)>3:
        if front_skipped: cross_points[1] = cross_points[3]
        if back_skipped: cross_points[-2] = cross_points[-4]

    for i in range(len(cross_points)-3,1,-2):
        if cross_points[i-1] == cross_points[i+1]:
            del cross_points[i:i+2]

    while len(cross_points)>2 and cross_points[2]<min_state_length:
        del cross_points[1:3]
    while len(cross_points)>2 and cross_points[-1]-cross_points[-3]<min_state_length:
        del cross_points[-3:-1]

    cross_points[-1] = chrom_length
    return cross_points


class _HMM(object):
    """A basic hidden markov model class for running the _HMM.gapped_viterbi method"""
    def __init__(self,states,priors,transition,observable,emission):
        self.states = states
        self.priors = [log(n) for n in priors]
        self.transition = [[log(n) for n in line] for line in transition]
        self.observable = {pair[1]:pair[0] for pair in enumerate(observable)}
        self.emission = [[log(n) for n in line] for line in emission]

    def get_most_likely_state(self,all_obs_points):
        obs_points = [snp for snp in all_obs_points if snp[1] in self.observable]
        enum_states = list(enumerate(self.states))
        state_sums = ((sum(self.emission[i][self.observable[snp[1]]] for snp in obs_points),state) for i,state in enum_states)
        best_state = max(state_sums,key=lambda x: x[0])[1]
        return best_state
        
    def gapped_viterbi(self,all_obs_points):
        vtrb = []
        trace_graph = {}

        obs_points = [snp for snp in all_obs_points if snp[1] in self.observable]
        if len(obs_points) < 1:
            warnings.warn("HMM defaulting to whole chromosome 'h' due to empty data.", UserWarning, stacklevel=2)
            return [0,"h",all_obs_points[-1][0]]
        
        enum_states = list(enumerate(self.states))

        for t,snp in enumerate(obs_points):
            obs = snp[1]
            vtrb.append([])
            if t==0:
                for i,state in enum_states:
                    vtrb[0].append(self.priors[i] + self.emission[i][self.observable[obs]])
                    trace_graph[(0,i)] = None
            else:
                gap_size = snp[0] - obs_points[t-1][0]
                for i,state in enum_states:
                    gap_transition_adjustment = log(gap_size)
                    max_prob,max_prev = max((vtrb[-2][p] + (self.transition[p][i]+gap_transition_adjustment) + self.emission[i][self.observable[obs]], p) for p,prev_state in enum_states)
                    vtrb[t].append(max_prob)
                    # 修复：这一行必须缩进在 for 循环内部！
                    trace_graph[(t,i)] = (t-1,max_prev)

        max_final,last = max((vtrb[-1][i],(len(vtrb)-1,i)) for i,state in enum_states)
        state_dets = []
        prev = last
        bp = []
        while prev!=None:
            state_dets.append(self.states[prev[1]])
            if trace_graph[prev] and trace_graph[prev][1]!=prev[1]:
                bp.append(prev[0])
            prev = trace_graph[prev]
        state_dets[:] = state_dets[::-1]
        bp[:] = bp[::-1]

        crosspoint_list = [0,state_dets[0]]
        for crosspoint in bp:
            crosspoint_list.append(obs_points[crosspoint][0])
            crosspoint_list.append(state_dets[crosspoint])
        crosspoint_list.append(all_obs_points[-1][0])

        return crosspoint_list

def _read_column(filename, col, filter=True):
    with open(filename, "r") as f:
        snp_list = []
        name = ""
        title_line = ""
        while title_line.strip()=="" or title_line.startswith("#"):
            title_line = f.readline().split("#")[0] 
        name = title_line.strip().split('\t')[col+2]
        for line in f:
            line = line.strip()
            if line!="" and not line.startswith("#"):
                line = line.split("#")[0] 
                items = line.split('\t')
                if filter==False or not items[col+2].lower()=="-": 
                    snp_list.append((int(float(items[1])), items[col+2].lower()))
        return snp_list,name

def _get_file_stats(filename):
    individual_count = 0
    with open(filename, "r") as f:
        title_line = ""
        while title_line.strip()=="" or title_line.startswith("#"):
            title_line = f.readline()
        indvs = title_line.strip().split('\t')[2:]
        individual_count = len(indvs)
    snps = _read_column(filename,0,filter=False)[0]
    last_index = snps[-1][0]
    snp_count = len(snps)
    return individual_count,snp_count,last_index
