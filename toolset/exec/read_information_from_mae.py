import glob,os

############# MODIFY HERE ############
## working directory names ##
dir_keywords = ['InducedFit_Lconstrain_enhanced_retain_no_sample25_AEKDEL_noH','InducedFit_Lconstrain_enhanced_retain_no_sample25_IGSLEL_noH','InducedFit_Lconstrain_enhanced_retain_no_sample25_IHSPDL_noH']
## key word for mae files ##
structure_keywords =['IFD_Lconstrain_enhanced_AEKDEL_noH_retain25_no_sample_ring','IFD_Lconstrain_enhanced_IGSLEL_noH_retain25_no_sample_ring','IFD_Lconstrain_enhanced_IHSPDL_noH_retain25_no_sample_ring']
#######################################

def get_information_from_mae(mae_file):
    f = open(mae_file,'r')
    lines=f.read().split('\n')
    f.close()
    tem=[]
    for i in range(len(lines)):
        if 'f_m_ct' in lines[i]:
            start_index= i
        if '}' in lines[i]:
            tem.append(i)
    f_m_ct_set = lines[start_index:tem[2]]
    tem=[]
    for i in range(len(f_m_ct_set)):
        if ':::' in f_m_ct_set[i]:
            tem.append(i)
        if '# First column is atom index #' in f_m_ct_set[i]:
            end_index = i-1

    name =f_m_ct_set[1:tem[0]]
    values =f_m_ct_set[tem[0]+1: end_index]
    dictionary={}                                                    #### create dictionary for storing energies
    for  i in range(len(name)):
      dictionary[name[i]] = values[i].split()[0]
    return dictionary

for j in range(len(dir_keywords)):
  f=open('IFDSCORE_'+dir_keywords[j]+'.txt','w')                     #### define output file name
  os.chdir(dir_keywords[j])
  FNs = glob.glob(structure_keywords[j]+'_*.mae')
  for i in range(1,len(FNs)+1):
    structure_name = structure_keywords[j]+'_'+str(i)+'.mae'
    information =get_information_from_mae(structure_name)
    if ' r_psp_IFDScore' in information.keys():                      #### extract IFD score
      IFD_score=information[' r_psp_IFDScore']
      print structure_name, IFD_score
      f.write(structure_name[:-4]+"%12.4f"%float(IFD_score)+'\n')
    else:
      print ('there is no IFD score for '+FNs[i])
  f.close()
  os.chdir('..')

